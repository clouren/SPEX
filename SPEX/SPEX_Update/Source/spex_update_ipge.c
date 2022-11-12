//------------------------------------------------------------------------------
// SPEX_Update/spex_update_ipge.c: perform one iteration of IPGE and perform
// any skipped scaling process.
//------------------------------------------------------------------------------

// TODO rewrite this function, try to extract common factor

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is to perform one iteration of IPGE. If the
// involving vector has pending scale factor(s), then it will be scaled and
// these factor(s) will be set to 1.
// In addition, this function should be used in successive IPGE
// updates, where history update could be involved for certain entry. Therefore,
// a history vector is required as a input/output. In case of a single IPGE
// iteration with no need for history update, this function should not be used
// (since the scaling for vector v can be skipped, refer to spex_update_dppu1).
// This function is called by the following functions:
// spex_update_dppu2: successive IPGE update for row k of U after swapping with
//                   row ks of U.
// spex_update_cppu: compute the (n-1)-th IPGE update of column k after vk is
//                   inserted.
// spex_update_triangular_solve: the REF triangular solve for LDx=v when L and
//                   v are sparse.
// spex_update_forward_sub: the forward substitution when solving LDUx=b, which
//                   is essensially REF triangular solve for LDx=b when b is
//                   dense.
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
// while x[perm[j+1:n]] need further update. All entries in vector x may have
// common factor x_scale that has not yet been applied to x[1:n].


#define SPEX_FREE_ALL                \
{                                    \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);            \
}

#include "spex_update_internal.h"

SPEX_info spex_update_ipge // perform IPGE on x based on v
(
    spex_scattered_vector sv_x, // array of size n for x in the scattered form.
                    // x could be dense if sv_x->i = NULL.
    int64_t *h,     // history vector for x, x[i] was last updated in the
                    // SPEX_FLIP(h[i])-th iteration
    int64_t *prev,  // prev is the index of the found previous entry of the last
                    // one (i.e., 2nd last entry) in v(perm). update if !prev
    SPEX_vector v,  // v is the vector that contains the j-th pivot
                    // used to compute x in the j-th IPGE iteration, which is
                    // the vector v in the equations mentioned above
    const int64_t *perm,        // permutation
    const int64_t *perm_inv,    // inverse of perm, can be NULL if prev == NULL
    const SPEX_matrix rhos,     // array of scaled pivots
    const int64_t j             // column index of v
)
{
    SPEX_info info;
    int sgn;
    // the first entry of v must be the pivot, which must be nonzero
    info = SPEX_mpz_sgn(&sgn, v->x[0]);
    if (info != SPEX_OK || sgn == 0 || v->i[0] != perm[j])
    {
        return SPEX_INCORRECT_INPUT;
    }

    mpq_t pending_scale; SPEX_MPQ_SET_NULL(pending_scale);
    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
    mpz_t *sd = rhos->x.mpz;
    int64_t perm_j = perm[j];

    // check if v has only 1 entry in the diagonal. If so, perform history
    // update if needed and return SPEX_OK
    int64_t real_hj = SPEX_FLIP(h[perm_j]);
    h[perm_j] = real_hj;
    if (v->nz == 1)
    {
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
        return SPEX_OK;
    }

    int64_t p, i, real_hi;
    int vscale = 0; // 1: v_scale == 1; -1: v_scale == -1; 0 none of above
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    SPEX_CHECK(SPEX_mpz_cmp_ui(&sgn, SPEX_MPQ_DEN(v->scale), 1));
    if (sgn == 0) // den(v_scale) == 1
    {
        SPEX_CHECK(SPEX_mpz_cmpabs_ui(&sgn, SPEX_MPQ_NUM(v->scale), 1));
        if (sgn == 0)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_MPQ_NUM(v->scale)));
            vscale = sgn < 0 ? -1 : 1;
        }
    }

    // pending_scale = x[perm[j]]/sd[h[perm[j]]]
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sv_x->x[perm_j]));
    if (real_hj > -1)
    {
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[real_hj]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
    }

    // The pivot in vector v is v->p[0], simple apply scale to it if needed
    if (vscale != 1)
    {
        SPEX_CHECK(SPEX_mpz_set(v->x[0], sd[j]));
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

        // apply scale factor to v(i): v(i) = v(i)*v_scale
        if (vscale == 0)
        {
            SPEX_CHECK(SPEX_mpz_divexact(v->x[p],
                                         v->x[p], SPEX_MPQ_DEN(v->scale)));
            SPEX_CHECK(SPEX_mpz_mul     (v->x[p],
                                         v->x[p], SPEX_MPQ_NUM(v->scale)));
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
            if (prev != NULL && perm_inv != NULL && perm_inv[i] > *prev &&
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
            SPEX_CHECK(SPEX_mpz_fdiv_q(tmpz, tmpz,
                                    SPEX_MPQ_DEN(pending_scale)));

            // x[i] = floor(x[i]/sd[h[i]])
            if (real_hi > -1)
            {
                SPEX_CHECK(SPEX_mpz_fdiv_q(sv_x->x[i], sv_x->x[i],sd[real_hi]));
            }

            // x[i] = x[i]- tmpz
            SPEX_CHECK(SPEX_mpz_sub(sv_x->x[i], sv_x->x[i], tmpz));
        }
        else
        {
            // -----------------------------------------------------------------
            // x[i] = (x[i]-v[i]*x[perm[j]])/sd[h[i]].
            // -----------------------------------------------------------------
            SPEX_CHECK(SPEX_mpz_submul(sv_x->x[i], v->x[p], sv_x->x[perm_j]));
            if (real_hi > -1)
            {
                SPEX_CHECK(SPEX_mpz_divexact(sv_x->x[i],
                                             sv_x->x[i], sd[real_hi]));
            }
        }

        // update h[i] and last_nz_b4_ks
        h[i] = SPEX_FLIP(j);
    }

    // reset v_scale to 1
    SPEX_CHECK(SPEX_mpq_set_ui(v->scale, 1, 1));

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
