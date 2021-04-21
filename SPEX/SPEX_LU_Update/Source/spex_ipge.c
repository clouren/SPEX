//------------------------------------------------------------------------------
//SPEX_CHOLMOD/spex_ipge.c: perform one iteration of IPGE and perform skipped
//                          any scaling process.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform one iteration of IPGE. If the
// involving vector has pending scale factor(s), then it will be scaled and
// these factor(s) will be set to 1.
// In addition, this function should be used in successive IPGE
// updates, where history update could be involved for certain entry. Therefore,
// a history vector is required as a input/output. In case of a single IPGE
// iteration with no need for history update, this function should not be used
// (since the scaling for vector v can be skipped, refer to spex_dppu1).
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
// Consider when performing the j-th IPGE update for x[i], and the j-th pivot
// is sd[j] and the j-th vector (row/column) is v. When all entries in v have
// no pending scale factor, v(perm[j])=sd[j], and we can perform IPGE for x[i]
// as
// x[i] = (x[i]*sd[j]-v[i]*x[perm[j]])/sd[j-1].
// This equation holds regardless of pending scaling factor x_scale is 1 or
// not.  However, when calling this function, x should have no pending scale.
// Otherwise, the result of IPGE update is not guaranteed to be in integer
// domain.
//
// In addition, in case of history update is needed before the IPGE update for
// x[i] and/or x[perm[j]], the equation becomes
// x[i] = x[i]*sd[j]/sd[h[i]]- v(i)*x[perm[j]]/sd[h[perm[j]]].
//
// When the IPGE update finished, all entries x[perm[1:j]] will be final,
// while x[perm[j+1:n]] need further update. All entries in vector x have common
// factor x_scale.


#define SPEX_FREE_ALL                \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_lu_update_internal.h"

SPEX_info spex_ipge // perform IPGE on x based on v
(
    spex_scattered_vector *sv_x,// array of size n for x in the scattered form.
                    // x could be dense if sv_x->i = NULL.
    int64_t *h,     // history vector for x, x[i] was last updated in the
                    // SPEX_FLIP(h[i])-th iteration
    int64_t *prev,  // prev is the index of the found previous entry of the last
                    // one (i.e., 2nd last entry) in v(perm). update if !prev
    SPEX_vector *v, // v is the vector that contains the j-th pivot
                    // used to compute x in the j-th IPGE iteration, which is
                    // the vector v in the equations mentioned above
    const int64_t *perm, // permutation
    const int64_t *perm_inv, // inverse of permutation
    const mpz_t *sd,// array of scaled pivots
    mpz_t *d,       // array of unscaled pivots
    const mpz_t new_dj,// new value for the j-th unscaled pivot
    mpq_t v_scale1, // the first pending scale for v
    mpq_t v_scale2, // the second pending scale for v
    mpq_t v_scale3, // a third pending scale not used for v
    const int64_t j, // column index of v
    const int64_t piv_j // the index of pivot in vector v//TODO remove
)
{
    SPEX_info info;
    if (!sv_x || !h || !perm || !perm_inv || !v || !sd || !d ||
        v->i[piv_j] != perm[j])
    {
        return SPEX_INCORRECT_INPUT;
    }

    mpq_t pending_scale; SPEX_MPQ_SET_NULL(pending_scale);// TODO make input
    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);

    // check if v has only 1 entry in the diagonal. If so, perform history
    // update if needed and return SPEX_OK
    int64_t real_hj = SPEX_FLIP(h[perm[j]]);
    h[perm[j]] = real_hj;
    if (v->nz == 1)
    {
        if (j-1 > real_hj) // require history update
        {
            SPEX_CHECK(SPEX_mpz_mul(sv_x->x[perm[j]],
                                    sv_x->x[perm[j]], sd[j-1]));
            if (real_hj > -1)
            {
                SPEX_CHECK(SPEX_mpz_divexact(sv_x->x[perm[j]],
                                             sv_x->x[perm[j]], sd[real_hj]));
            }
        }
        return SPEX_OK;
    }

    int64_t p, i, real_hi;
    int sgn;
    int vscale = 0; // 1: v_scale1 == 1; -1: v_scale1 == -1; 0 none of above
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    // v_scale3 = v_scale3*v_scale2
    SPEX_CHECK(SPEX_mpq_mul(v_scale3, v_scale3, v_scale2));
    // v_scale1 = v_scale1*v_scale2
    SPEX_CHECK(SPEX_mpq_mul(v_scale1, v_scale1, v_scale2));
    SPEX_CHECK(SPEX_mpq_set_ui(v_scale2, 1, 1));
    SPEX_CHECK(SPEX_mpz_cmp_ui(&sgn, SPEX_MPQ_DEN(v_scale1), 1));
    if (sgn == 0) // den(v_scale1) == 1
    {
        SPEX_CHECK(SPEX_mpz_cmpabs_ui(&sgn, SPEX_MPQ_NUM(v_scale1), 1));
        if (sgn == 0)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_MPQ_NUM(v_scale1)));
            vscale = sgn < 0 ? -1 : 1;
        }
    }

    // pending_scale = x[perm[j]]/sd[h[perm[j]]]
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sv_x->x[perm[j]]));
    if (real_hj > -1)
    {
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[real_hj]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
    }

    // NOTE: this could cause fillin in x 
    for (p = 0; p < v->nz; p++)
    {
        // no need to update x(perm[j]) but apply scale to v[perm[j]] if needed
        if (p == piv_j) // same as (i == perm[j])
        {
            if (vscale != 1)
            {
                SPEX_CHECK(SPEX_mpz_set(v->x[p], sd[j]));
                SPEX_CHECK(SPEX_mpz_set(d[j], new_dj));
            }
            continue;
        }
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[p]));
        if (sgn == 0)    // v[i] == 0
        {
            continue;
        }
        // column/row index in v
        i = v->i[p];
        real_hi = SPEX_FLIP(h[i]);

        // apply scale factor to v(i): v(i) = v(i)*v_scale1
        if (vscale == 0)
        {
            SPEX_CHECK(SPEX_mpz_divexact(v->x[p],
                                         v->x[p], SPEX_MPQ_DEN(v_scale1)));
            SPEX_CHECK(SPEX_mpz_mul     (v->x[p],
                                         v->x[p], SPEX_MPQ_NUM(v_scale1)));
        }
        else if (vscale == -1)
        {
            SPEX_CHECK(SPEX_mpz_neg(v->x[p], v->x[p]));
        }

        // x[i] = x[i]*sd[j]
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[i]));
        if (sgn != 0)    // x[i] != 0
        {
            // x[i] = x[i]*sd[j]
            SPEX_CHECK(SPEX_mpz_mul(sv_x->x[i], sv_x->x[i], sd[j]));
        }
        else if (sv_x->i != NULL && h[i] >= -1)
        {
            ASSERT(sv_x->nz < sv_x->nzmax);
            // this entry was not in nnz pattern, so we add it to nnz pattern
            sv_x->i[sv_x->nz] = i;
            sv_x->nz ++;

            // update prev if needed
            if (prev != NULL && perm_inv[i] > *prev &&
                perm_inv[i] != sv_x->nzmax-1)
            {
                *prev = perm_inv[i];
            }
        }
        if (real_hi != real_hj)
        {
            // -----------------------------------------------------------------
            // x[i] = x[i]/sd[h[i]]- v(i)*x[perm[j]]/sd[h[perm[j]]].
            // -----------------------------------------------------------------
            // tmpz = floor(v(i)*pending_scale)
            SPEX_CHECK(SPEX_mpz_mul(tmpz, v->x[p],
                                    SPEX_MPQ_NUM(pending_scale)));
#ifndef SPEX_DEBUG
            SPEX_CHECK(SPEX_mpz_fdiv_q(tmpz, tmpz,
                                    SPEX_MPQ_DEN(pending_scale)));

            // x[i] = floor(x[i]/sd[h[i]])
            if (real_hi > -1)
            {
                SPEX_CHECK(SPEX_mpz_fdiv_q(sv_x->x[i], sv_x->x[i],sd[real_hi]));
            }
#else
                mpq_t tmpq1; mpq_init(tmpq1);
                mpq_t tmpq2; mpq_init(tmpq2);mpq_set_ui(tmpq2,0,1);
            mpz_fdiv_qr(tmpz, SPEX_MPQ_NUM(tmpq1), tmpz, SPEX_MPQ_DEN(pending_scale));
                mpq_set_den(tmpq1,SPEX_MPQ_DEN(pending_scale));
                mpq_canonicalize(tmpq1);

            // x[i] = floor(x[i]/sd[h[i]])
            if (real_hi > -1)
            {
                mpz_fdiv_qr(sv_x->x[i], SPEX_MPQ_NUM(tmpq2), sv_x->x[i],sd[real_hi]);
                mpq_set_den(tmpq2,sd[real_hi]);
                mpq_canonicalize(tmpq2);
            }

            SPEX_CHECK(SPEX_mpq_cmp(&sgn,tmpq1,tmpq2));
            if (sgn!=0)
            {
                printf("file %s line %d\n",__FILE__,__LINE__);
            mpq_clear(tmpq1);
            mpq_clear(tmpq2);
                SPEX_CHECK(SPEX_PANIC);
            }
            mpq_clear(tmpq1);
            mpq_clear(tmpq2);
#endif

            // x[i] = x[i]- tmpz
            SPEX_CHECK(SPEX_mpz_sub(sv_x->x[i], sv_x->x[i], tmpz));
        }
        else
        {
            // -----------------------------------------------------------------
            // x[i] = (x[i]-v[i]*x[perm[j]])/sd[h[i]].
            // -----------------------------------------------------------------
            mpz_t xi; mpz_init_set(xi,sv_x->x[i]);
            SPEX_CHECK(SPEX_mpz_submul(sv_x->x[i], v->x[p], sv_x->x[perm[j]]));
            if (real_hi > -1)
            {
                SPEX_CHECK(SPEX_mpz_divexact(sv_x->x[i],
                                             sv_x->x[i], sd[real_hi]));
            }
        }

        // update h[i] and last_nz_b4_ks
        h[i] = SPEX_FLIP(j);
    }

    // reset v_scale1 to 1
    SPEX_CHECK(SPEX_mpq_set_ui(v_scale1, 1, 1));

    if (j-1 > real_hj) // require history update
    {
        SPEX_CHECK(SPEX_mpz_mul(sv_x->x[perm[j]],
                                sv_x->x[perm[j]], sd[j-1]));
        if (real_hj > -1)
        {
            SPEX_CHECK(SPEX_mpz_divexact(sv_x->x[perm[j]],
                                         sv_x->x[perm[j]], sd[real_hj]));
        }
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
