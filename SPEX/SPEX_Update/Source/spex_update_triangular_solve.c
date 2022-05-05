//------------------------------------------------------------------------------
// SPEX_Update/spex_update_triangular_solve.c: perform REF triangular solve up
// to specified iteration, additional history update for certain entries should
// be done after calling this function.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is to perform REF triangular solve for LDx=v up to
// specified IPGE iteration when both L and v are sparse. Additional history
// update should be done for certain entries based on the history vector. This
// function is used internally when computing the inserted column for U. For
// forward solving LDUx=b where b is mostly considered and stored as dense
// vector, we will use spex_forward_sub.

#include "spex_update_internal.h"

SPEX_info spex_update_triangular_solve // perform REF triangular solve for LDx=v
(
    spex_scattered_vector *sv_x,// the scattered version of solution for LDx=v,
                        // using the first k-1 columns of L
    int64_t *x_top,     // P_inv[sv_x->i[0...(*x_top)]] <= (*last_update), that
                        // is, sv_x->i[0...(*x_top)] give the indices of all
                        // entries that are up-to-date. However, this is updated
                        // only when i_2ndlast is requested.
    int64_t *h,         // history vector for x
    int64_t *last_update,// the number of finished IPGE iterations, which is
                        // also the number of columns in L used last time
    int64_t *i_2ndlast, // i_2ndlast is the index of the found last nnz entry
                        // of x[P] in the range of (last_update, n-2], this
                        // could be NULL if not needed
    const int64_t k,    // compute x up to k-th IPGE iteration, that is, using
                        // the first k-1 columns of L
    const SPEX_matrix *L,  // matrix L
    const SPEX_matrix *U,  // matrix U
    const SPEX_matrix *rhos,// array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
)
{
    SPEX_info info;
    int sgn;
    int64_t j;
    int64_t n = sv_x->nzmax;

    // there is no nnz in vk(P[last_update,n-2])
    if (i_2ndlast != NULL && *i_2ndlast == -1)
    {
        return SPEX_OK;    
    }

    if (*last_update < k-1)
    {
        // iterate across each entry
        for (j = *last_update+1; j < k; j++)
        {
            // skip if x(P[j]) == 0
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[P[j]]));
            if (sgn == 0)       {continue; }

            // perform j-th IPGE update for x
            SPEX_CHECK(spex_update_ipge(sv_x, h, i_2ndlast, L->v[j], P,
                P_inv, rhos, j));
        }
        *last_update = k-1;

        // double check if i_2ndlast gives the correct index for 2nd last nnz
        if(i_2ndlast != NULL && *i_2ndlast != -1)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[*i_2ndlast]));
            if (sgn == 0 || *i_2ndlast < k) // needs to update i_2ndlast
            {
                int64_t p, real_j;
                *i_2ndlast = -1;
                for (p = *x_top; p < sv_x->nz;)
                {
                    j = sv_x->i[p];
                    real_j = P_inv[j];

                    // skip entries above k-th row
                    if (real_j <= *last_update)
                    {
                        // move these entries before x_top
                        sv_x->i[p] = sv_x->i[*x_top];
                        sv_x->i[*x_top] = j;
                        (*x_top)++;
                        p++;
                        continue;
                    }

                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[j]));
                    if (sgn == 0)
                    {
                        // remove it from nnz pattern
                        h[j] = -1;
                        sv_x->nz--;
                        sv_x->i[p] = sv_x->i[sv_x->nz];
                    }
                    else
                    {
                        if (real_j > *i_2ndlast && real_j != n-1)
                        {
                            *i_2ndlast = real_j;
                        }
                        p++;
                    }
                }
            }
        }
    }

    return SPEX_OK;
}
