//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_triangular_solve.c: perform REF triangular solve up to
// specified iteration, additional history update for certain entries should
// be done after calling this function.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform REF triangular solve for LDx=v up to
// specified IPGE iteration when both L and v are sparse. Additional history
// update should be done for certain entries based on the history vector. This
// function is used internally when computing the inserted column for U. For
// forward solving LDUx=b where b is mostly considered and stored as dense
// vector, we will use spex_forward_sub.

#include "spex_lu_update_internal.h"

SPEX_info spex_triangular_solve // perform REF triangular solve for LDx=v
(
    spex_scattered_vector *sv_x,// the scattered version of solution for LDx=v,
                        // using the first k-1 columns of L
    int64_t *h,         // history vector for x
    int64_t *last_update,// the number of finished IPGE iterations, which is
                        // also the number of columns in L used last time
    int64_t *i_2ndlast, // i_2ndlast is the index of the found last nnz entry
                        // of x[P] in the range of (last_update, n-2], this
                        // could be NULL if not needed
    const int64_t k,    // compute x up to k-th IPGE iteration, that is, using
                        // the first k-1 columns of L
    const SPEX_mat *L,// matrix L
    const SPEX_mat *U,// matrix U
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Ucp, // col pointers for col-wise nnz pattern of U
    const int64_t *Ucx, // the value of k-th entry is found as
                        // U->v[Uci[k]]->x[Ucx[k]]
    mpq_t *S,           // the pending scale factor matrix
    const mpz_t *sd,    // array of scaled pivots
    mpz_t *d,           // array of unscaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t *Q    // column permutation
)
{
    SPEX_info info;
    int sgn;
    int64_t j;
    if (!sv_x || !h || !last_update || !L || !S || !P || !P_inv || !sd)
    {
        return SPEX_INCORRECT_INPUT;
    }

    // there is no nnz in vk(P[last_update,n-2]), so proceed only when k=n-1
    int64_t n = sv_x->nzmax;
    if (i_2ndlast != NULL && *i_2ndlast == -1 && k != n-1)
    {
        return SPEX_OK;    
    }

    if (*last_update < k-1)
    {
        // TODO iterate sorted nnz pattern?
        for (j = *last_update+1; j < k; j++)
        {
            // skip if x(P[j]) == 0
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[P[j]]));
            if (sgn == 0)       {continue; }

            // TODO add this to all caller of spex_ipge,
            // but no need to update spex_ipge
            ASSERT(L->v[j]->i[Ldiag[j]] == P[j]);
            // perform j-th IPGE update for x
            SPEX_CHECK(spex_ipge(sv_x, h, i_2ndlast, L->v[j], P, P_inv,
                sd, d, U->v[j]->x[Ucx[Ucp[Q[j]+1]-1]],
                SPEX_LUU_2D(S, 1, j), SPEX_LUU_2D(S, 3, j), SPEX_LUU_2D(S, 2, j),
                j, Ldiag[j]));
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
                for (p = 0; p < sv_x->nz;)
                {
                    j = sv_x->i[p];
                    real_j = P_inv[j];

                    // skip entries above k-th row
                    if (real_j <= *last_update)
                    {
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
