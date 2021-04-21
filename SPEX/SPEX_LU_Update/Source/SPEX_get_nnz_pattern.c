//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_get_nnz_pattern.c: get the row-wise nonzero pattern of L
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
    SPEX_FREE(Ldiag_new);                  \
    SPEX_FREE(Lr_offdiag_new);             \
    SPEX_FREE(Ucp_new);                    \
    SPEX_FREE(Uci_new);                    \
    SPEX_FREE(Ucx_new);

#include "spex_lu_update_internal.h"

SPEX_info SPEX_get_nnz_pattern    // find the nnz pattern of L and U
(
    // OUTPUT:
    int64_t **Ldiag,              // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    int64_t **Lr_offdiag,         // Lr_offdiag[k] gives the column index of the
                                  // last off-diagonal nnz in k-th row of L.
                                  // -1 if no off diagonal entry
    int64_t **Uci,                // the row index for col-wise nnz pattern of U
    int64_t **Ucp,                // col pointers for col-wise pattern of U
    int64_t **Ucx,                // find the value of k-th entry as
                                  // U->v[Uci[k]]->x[Ucx[k]]
    // INPUT:
    const SPEX_mat *L,         // the target matrix L
    const SPEX_mat *U,         // the target matrix U
    const int64_t *P,             // row permutation
    const SPEX_options *option    // command option
)
{
    // inputs are checked in SPEX_LUU
    SPEX_info info;
    int64_t *Ldiag_new = NULL, *Lr_offdiag_new = NULL, *Ucp_new = NULL,
            *Uci_new = NULL, *Ucx_new = NULL, *rowcount = NULL;
    int64_t i, j, p;
    int64_t n = U->n;

    rowcount       = (int64_t*) SPEX_calloc( n,   sizeof(int64_t));
    Ldiag_new      = (int64_t*) SPEX_malloc( n   *sizeof(int64_t));
    Lr_offdiag_new = (int64_t*) SPEX_malloc( n   *sizeof(int64_t));
    Ucp_new        = (int64_t*) SPEX_malloc((n+1)*sizeof(int64_t));
    if (!rowcount || !Ldiag_new || !Lr_offdiag_new || !Ucp_new)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    for (i = 0 ; i < n ; i++)
    {
        Lr_offdiag_new[i] = -1;
        Ldiag_new[i] = -1;
    }
    for (j = 0; j < n; j++)
    {
        for (p = 0; p < L->v[j]->nz; p++)
        {
            // row index
            i = L->v[j]->i[p];

            if (i == P[j])  // current entry is diagonal of L(P,:)
            {
#ifdef SPEX_DEBUG
                int sgn;
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[j]->x[p]));
                ASSERT (sgn != 0);
#endif
                Ldiag_new[j] = p; // get the row pointer
            }
            else            // get the last row-wise off-diagonal entries
            {
                Lr_offdiag_new[i] = j; // get the column index
            }
        }
        if (Ldiag_new[j] == -1)
        {
            // column j of L has no nnz
            return SPEX_INCORRECT_INPUT;
        }
    }

    for (i = 0; i < n; i++)
    {
        for (p = 0 ; p < U->v[i]->nz; p++)
        {
            j = U->v[i]->i[p];
            rowcount [j]++ ;
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

    for (i = 0 ; i < n ; i++)
    {
        for (p = 0 ; p < U->v[i]->nz ; p++)
        {
            j = U->v[i]->i[p];
            Uci_new[ rowcount[j] ] = i ;
            Ucx_new[ rowcount[j] ] = p ;
            rowcount[j] ++;
        }
    }

    *Ldiag = Ldiag_new;
    *Lr_offdiag = Lr_offdiag_new;
    *Ucp = Ucp_new;
    *Uci = Uci_new;
    *Ucx = Ucx_new;
    SPEX_FREE_WORK;
    return SPEX_OK;
}
