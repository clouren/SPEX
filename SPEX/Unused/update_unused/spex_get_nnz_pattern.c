//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_get_nnz_pattern.c: get the row-wise nonzero pattern of L
// and column-wise nonzero pattern of U.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function finds the column-wise nonzero pattern from a
// compressed-row matrix U, and the column index of last off-diagonal entry in
// each row and the pointer to each diagonal entry of a compressed-column
// matrix L. 

#define SPEX_FREE_WORK                     \
    SPEX_FREE(rowcount);

#define SPEX_FREE_ALL                      \
    SPEX_FREE_WORK;                        \
    SPEX_FREE(Lr_offdiag_new);             \
    SPEX_FREE(Uc_offdiag_new);             \
    SPEX_FREE(Ucp_new);                    \
    SPEX_FREE(Uci_new);                    \
    SPEX_FREE(Ucx_new);

#include "spex_lu_update_internal.h"

SPEX_info spex_get_nnz_pattern    // find the nnz pattern of L and U
(
    // OUTPUT:
    int64_t **Lr_offdiag,         // Lr_offdiag[k] gives the column index of the
                                  // last off-diagonal nnz in k-th row of L.
                                  // -1 if no off diagonal entry
    int64_t **Uc_offdiag,         // Lr_offdiag[k] gives the column index of the
    int64_t **Uci,                // the row index for col-wise nnz pattern of U
    int64_t **Ucp,                // col pointers for col-wise pattern of U
    int64_t **Ucx,                // find the value of k-th entry as
                                  // U->v[Uci[k]]->x[Ucx[k]]
    // INPUT:
    const SPEX_mat *L,            // the target matrix L
    const SPEX_mat *U,            // the target matrix U
    const int64_t *P_inv,         // inverse of row permutation
    const int64_t *Q_inv,         // inverse of col permutation
    const int64_t k,              // index of column to be replaced
    const SPEX_options *option    // command option
)
{
    // inputs are checked in SPEX_LUU
    SPEX_info info;
    int64_t *Lr_offdiag_new = NULL, *Uc_offdiag_new = NULL, *Ucp_new = NULL,
            *Uci_new = NULL, *Ucx_new = NULL, *rowcount = NULL;
    int64_t i, j, p;
    int64_t n = U->n;

    rowcount       = (int64_t*) SPEX_calloc( n,   sizeof(int64_t));
    Lr_offdiag_new = (int64_t*) SPEX_malloc( n   *sizeof(int64_t));
    Ucp_new        = (int64_t*) SPEX_malloc((n+1)*sizeof(int64_t));
    Uc_offdiag_new = (int64_t*) SPEX_malloc( n   *sizeof(int64_t));
    if (!rowcount || !Lr_offdiag_new || !Ucp_new)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    int64_t Q_inv_k = Q_inv[k];
    // initialize entries Q_inv_k+1:n
    for (i = Q_inv_k+1; i < n ; i++)
    {
        Lr_offdiag_new[i] = -1;
        Uc_offdiag_new[i] = -1;
    }
    for (j = Q_inv_k+1; j < n; j++)
    {
#ifdef SPEX_DEBUG
        // the first entry should be the corresponding pivot
        ASSERT(P_inv[L->v[j]->i[0]] == j);
        int sgn;
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[j]->x[0]));
        ASSERT (sgn != 0);
#endif
        for (p = 1; p < L->v[j]->nz; p++)
        {
            // row index
            i = P_inv[L->v[j]->i[p]];

            // get the last row-wise off-diagonal entries
            Lr_offdiag_new[i] = j; // get the column index
        }
    }

    // remove column k of U
    for (i = 0; i < Q_inv_k; i++)
    {
        for (p = 1 ; p < U->v[i]->nz; p++)
        {
            j = U->v[i]->i[p];
            if (j == k)
            {
                // move the last entry to current position
                U->v[i]->nz--;
                SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[p],
                                         U->v[i]->x[U->v[i]->nz]));
                U->v[i]->i[p] = U->v[i]->i[U->v[i]->nz];
                break;
            }
        }
    }
    for (i = Q_inv_k+1; i < n; i++)
    {
        ASSERT(Q_inv[U->v[i]->i[0]] == i);
        for (p = 0 ; p < U->v[i]->nz; p++)
        {
            j = U->v[i]->i[p];
            rowcount [j]++ ;
            if (p != 0)
            {
                Uc_offdiag_new[Q_inv[j]] = i;
            }
        }
    }

    // compute cumulative sum of rowcount to get the col pointer
    SPEX_CHECK(SPEX_cumsum(Ucp_new, rowcount, n));

    int64_t U_nnz = Ucp_new[n];
    Uci_new = (int64_t*) SPEX_malloc(U_nnz*sizeof(int64_t));
    Ucx_new = (int64_t*) SPEX_malloc(U_nnz*sizeof(int64_t));
    if (!Uci_new || !Ucx_new)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    for (i = Q_inv_k+1 ; i < n ; i++)
    {
        for (p = 0 ; p < U->v[i]->nz ; p++)
        {
            j = U->v[i]->i[p];
            Uci_new[ rowcount[j] ] = i ;
            Ucx_new[ rowcount[j] ] = p ;
            rowcount[j] ++;
        }
    }

    *Lr_offdiag = Lr_offdiag_new;
    *Uc_offdiag = Uc_offdiag_new;
    *Ucp = Ucp_new;
    *Uci = Uci_new;
    *Ucx = Ucx_new;
    SPEX_FREE_WORK;
    return SPEX_OK;
}
