//------------------------------------------------------------------------------
//SPEX_CHOLMOD/spex_ipge.c: perform one iteration of IPGE with effort to skip
//                          possible scaling process.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform one iteration of IPGE when the
// involving vectors have pending scale factors that are wished not to be
// applied. In addition, this function should be used in successive IPGE
// update, where history update could be involved for certain entry. Therefore,
// a history vector is required as a input/output. In case of a single IPGE
// iteration with no need for history update, this function should not be used.
// This function is called by the following functions:
// spex_dppu2: successive IPGE update for row k of U after swapping with row
//             ks of U.
// spex_cppu: compute the (n-1)-th IPGE update of column k after vk is inserted.
// spex_triangular_solve: the REF triangular solve for LDx=v when L and v are
//             sparse.
// spex_forward_sub: the forward substitution when solving LDUx=b, which is
//             essensially REF triangular solve for LDx=b when b is dense.
//
// Algorithm explanation:
// The input matrix M can be L or U, and perm should be correspondingly P or Q.
// To make the explaination simpler, we assume M = L, perm = P and
// perm_inv = P_inv. The case for U can be easily generalized.
//
// Consider when performing the j-th IPGE update for x[i], and the j-th col of
// L is L(:,j) = L->v[j] and its pivot is L(P[j],j) = L->v[j]->x[0]. When all
// entries in L->v[j] have no pending scale factor, i.e., S(1,j) = 1,
// L(P[j],j) = L->v[j]->x[0] = sd[j],
// and we can perform IPGE for x[i] as
// x[i] = (x[i]*sd[j]-L(i,j)*x[P[j]])/sd[j-1].
//
// However, in the case when S(1,j) = L->v[j]->scale != 1, and we want to keep
// it as it is, we can compute x[i] as
// x[i] = (x[i]*L(P[j],j)-L(i,j)*x[P[j]])/L(P[j-1],j-1)
// since L(P[j],j)*S(1,j)=sd[j]. When the update finished, entries that are
// updated in the j-th IPGE (i.e., x[P[j+1, ..., n-1]]) need to multiply with
// S(1,j)/S(1,j-1), while x_scale is maintained unchanged in this function.
// However, this should be handled by the caller of this function.
// Simplification can be done to compute the final value of each entry of x,
// instead of having them updated at the end of each call of this function.
// An efficient way of getting the final value is provided below.
//
// In addition, in case of history update is needed before the IPGE update for
// x[i] and/or x[P[j]], the equation becomes
// x[i] = x[i]*L(P[j],j)/L(P[h[i]],h[i])- L(i,j)*x[P[j]]/L(P[h[P[j]]], h[P[j]]).
//
// As mentioned above, this function is called to perform successive IPGE
// updates. Let us denote these IPGE updates as [js, js+1, ..., je-1, je]-th
// updates, then the exact final value of each entry of x can be computed as
//  # x[P[0, 1, ..., js]] should be obtained by multiplying with x_scale, which
//    is NOT changed after each call of this function, which is
//    x[P[0, 1, ..., js]]   *= x_scale
//  # x[P[js+1]]            *= x_scale*S(1,js)  /S(1,js-1)
//  # x[P[js+2]]            *= x_scale*S(1,js+1)/S(1,js-1)
//    ...
//  # x[P[je]]              *= x_scale*S(1,je-1)/S(1,js-1)
//  # x[P[je+1, ..., n-1]]  *= x_scale*S(1,je)  /S(1,js-1)
// Therefore, the caller of this function should make sure all entries share
// the same pending scale.


#define SPEX_FREE_ALL                \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_lu_update_internal.h"

SPEX_info spex_ipge1 // perform IPGE on x based on v
(
    spex_scattered_vector *sv_x,// array of size n for x in the scattered form.
                    // x could be dense by setting sv_x->i = NULL.
    int64_t *h,     // history vector for x, x[i] was last updated in the
                    // SPEX_FLIP(h[i])-th iteration
    int64_t *prev,  // prev is the index of the found previous entry of the last
                    // one (i.e., 2nd last entry) in v(perm). update if !prev
    const SPEX_mat *M,// M->v[j] is the vector that contains the j-th pivot
                    // used to compute x in the j-th IPGE iteration, which is
                    // the vector v in the equations mentioned above
    const int64_t *perm, // permutation
    const int64_t *perm_inv, // inverse of permutation
    const mpz_t *sd,
    mpq_t v_scale,
    const int64_t j // column index of v in M
)
{
    SPEX_info info;
    if (!sv_x || !h || !perm || (prev != NULL && !perm_inv) || !M)
    {
        return SPEX_INCORRECT_INPUT;
    }

    int64_t p, i, real_hi, perm_j = perm[j];
    int sgn;
    int64_t real_hj = SPEX_FLIP(h[perm_j]);
    h[perm_j] = real_hj;

    SPEX_vector *v = M->v[j];
    mpq_t pending_scale; SPEX_MPQ_SET_NULL(pending_scale);
    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    // pending_scale = x[perm[j]]/M->v[h[perm[j]]]->x[0]
    SPEX_CHECK(SPEX_mpq_set(pending_scale, v_scale));
    if (real_hj > -1)
    {
        SPEX_CHECK(SPEX_mpz_mul(SPEX_MPQ_DEN(pending_scale),
                                SPEX_MPQ_DEN(pending_scale), sd[real_hj]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
    }
    // perform IPGE for x if the corresponding entry in v != 0, skip x(perm[j])
    // NOTE: this could cause fillin in x 
    for (p = 1; p < v->nz; p++)
    {
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[p]));
        if (sgn == 0)    // v[i] == 0
        {
            continue;
        }
        // column/row index in v
        i = v->i[p];
        real_hi = SPEX_FLIP(h[i]);

        SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[i]));
        if (sgn != 0)    // x[i] != 0
        {
            if (real_hi == real_hj)
            {
                // x[i] = x[i]*v[perm[j]]
                SPEX_CHECK(SPEX_mpz_mul(sv_x->x[i], sv_x->x[i], v->x[0]));
            }
            else
            {
                // x[i] = x[i]*sd[j]
                SPEX_CHECK(SPEX_mpz_mul(sv_x->x[i], sv_x->x[i], sd[j]));
            }
        }
        else if (sv_x->i != NULL && h[i] >= -1)
        {
            // this entry was not in nnz pattern, so we add it to nnz pattern
            sv_x->i[sv_x->nz] = i;
            sv_x->nz ++;

            // update prev if needed
            if (prev != NULL && perm_inv != NULL && perm_inv[i] > *prev &&
                perm_inv[i] != sv_x->nzmax-1)
            {
                *prev = perm_inv[i];
            }
        }
        if (real_hi != real_hj)
        {
            // -----------------------------------------------------------------
            // x[i] = x[i]/M->v[h[i]]->x[0]- v(i)*pending_scale
            // -----------------------------------------------------------------
//      = x[i]*sd[j]/sd[h[i]]             -L(i,j)*S(1,j)*x[P[j]]/sd[h[P[j]]]
            // tmpz = floor(v(i)*pending_scale)
            SPEX_CHECK(SPEX_mpz_mul(tmpz, v->x[p], sv_x->x[perm_j]));
            SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz,
                                    SPEX_MPQ_NUM(pending_scale)));
#ifndef SPEX_DEBUG
            SPEX_CHECK(SPEX_mpz_fdiv_q(tmpz, tmpz,
                                    SPEX_MPQ_DEN(pending_scale)));

            // x[i] = x[i]/sd[h[i]]
            if (real_hi > -1)//TODO move into check if x[i]==0
            {
                SPEX_CHECK(SPEX_mpz_fdiv_q(sv_x->x[i],
                                           sv_x->x[i], sd[real_hi]));
            }
#else
                mpq_t tmpq1; mpq_init(tmpq1);
                mpq_t tmpq2; mpq_init(tmpq2);mpq_set_ui(tmpq2,0,1);
                mpz_fdiv_qr(tmpz, SPEX_MPQ_NUM(tmpq1),
                                            tmpz, SPEX_MPQ_DEN(pending_scale));
                SPEX_CHECK(SPEX_mpq_set_den(tmpq1,SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpq_canonicalize(tmpq1));

            // x[i] = floor(x[i]/sd[h[i]])
            if (real_hi > -1)
            {
                mpz_fdiv_qr(sv_x->x[i], SPEX_MPQ_NUM(tmpq2),
                            sv_x->x[i], sd[real_hi]);
                SPEX_CHECK(SPEX_mpq_set_den(tmpq2,sd[real_hi]));
                SPEX_CHECK(SPEX_mpq_canonicalize(tmpq2));
            }

            SPEX_CHECK(SPEX_mpq_cmp(&sgn,tmpq1,tmpq2));
            if (sgn!=0)
            {
                printf("%ld-th IPGE at p=%ld, perm_inv[i] = perm_inv[%ld]=%f\n",
                    j,p,i,(perm_inv == NULL)?NAN:(double)perm_inv[i]);
                    gmp_printf("%Zd\n%Zd\n",SPEX_MPQ_NUM(tmpq1),SPEX_MPQ_DEN(tmpq1));
                SPEX_CHECK(SPEX_mpz_divexact(tmpz,
                                             v->x[p], SPEX_MPQ_DEN(v_scale)));
                SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(v_scale)));
                SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz, sv_x->x[perm_j]));
                if (real_hj > -1)
                {
                    mpz_fdiv_qr(tmpz, SPEX_MPQ_NUM(tmpq1),
                                tmpz, sd[real_hj]);
                    SPEX_CHECK(SPEX_mpq_set_den(tmpq1, sd[real_hj]));
                    SPEX_CHECK(SPEX_mpq_canonicalize(tmpq1));
                }
                else
                {
                    gmp_printf("%Qd\n",tmpq1);
                    mpq_set_ui(tmpq1, 0,1);
                }
                    gmp_printf("%Qd\n",tmpq1);
            SPEX_CHECK(SPEX_mpq_cmp(&sgn,tmpq1,tmpq2));
            printf("%s\n", sgn == 0?"same":"different");
                
                SPEX_MPQ_CLEAR(tmpq1);
                SPEX_MPQ_CLEAR(tmpq2);
                SPEX_CHECK(SPEX_PANIC);
            }
            SPEX_MPQ_CLEAR(tmpq1);
            SPEX_MPQ_CLEAR(tmpq2);

#endif

            // x[i] = x[i]- tmpz
            SPEX_CHECK(SPEX_mpz_sub(sv_x->x[i], sv_x->x[i], tmpz));
        }
        else
        {
            // -----------------------------------------------------------------
            // x[i] = (x[i]-v[i]*x[perm[j]])/M->v[h[i]]->x[0].
            // -----------------------------------------------------------------
//      = (x[i]*L(P[j],j)-L(i,j),x[P[j]])*S(1,j)/sd[h[i]]
            SPEX_CHECK(SPEX_mpz_submul(sv_x->x[i], v->x[p], sv_x->x[perm_j]));
            SPEX_CHECK(SPEX_mpz_divexact(sv_x->x[i],
                                         sv_x->x[i],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sv_x->x[i], sv_x->x[i],
                                    SPEX_MPQ_NUM(pending_scale)));
        }

        // update h[i] and last_nz_b4_ks
        h[i] = SPEX_FLIP(j);
    }

    if (j-1 > real_hj) // require history update
    {
        SPEX_CHECK(SPEX_mpz_mul(sv_x->x[perm_j],
                                sv_x->x[perm_j], sd[j-1]));
        if (real_hj > -1)
        {
            SPEX_CHECK(SPEX_mpz_divexact(sv_x->x[perm_j],
                                         sv_x->x[perm_j], sd[real_hj]));
        }

    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
// x[i] = (x[i]*L(P[j],j)-L(i,j)*x[P[j]])*S(1,j)/sd(j-1)

// x[i] = (x[i]*sd[j-1]/sd[h[i]]*L(P[j],j)-L(i,j)*x[P[j]])*S(1,j)/sd[j-1]
//      = x[i]*sd[j-1]/sd[h[i]]*L(P[j],j)*S(i,j)/sd[j-1]
//                                        -L(i,j)*x[P[j]]*S(1,j)/sd[j-1]
//      = x[i]*sd[j]/sd[h[i]]             -L(i,j)*S(1,j)*x[P[j]]/sd[h[P[j]]]
//      = (x[i]*L(P[j],j)-L(i,j),x[P[j]])*S(1,j)/sd[h[i]]

// x[i] = x[i]*L(P[j],j)/L(P[h[i]],h[i])- L(i,j)*x[P[j]]/L(P[h[P[j]]], h[P[j]]).
