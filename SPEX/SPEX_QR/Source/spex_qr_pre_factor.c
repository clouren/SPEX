//------------------------------------------------------------------------------
// SPEX_QR/spex_qr_pre_factor: Symbolic left-looking Cholesky for QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Lorena Mejia Domenzain, Christopher Lourenco,
// Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_FREE(w);                   \
}

# define SPEX_FREE_ALL               \
{                                    \
    SPEX_FREE_WORKSPACE              \
    SPEX_matrix_free(&R, NULL);      \
}

#include "spex_cholesky_internal.h"


/* Purpose: This function performs a symbolic sparse triangular solve for
 * each column of R
 * It allocates the memory for the R matrix and determines the full nonzero
 * pattern of R
 *
 * Importantly, this function assumes that A has already been permuted.
 *
 * Input arguments of the function:
 *
 * R_handle:    A handle to the R matrix. Null on input.
 *              On output, contains a pointer to the partial R matrix.
 *
 * A:           The user's permuted input matrix
 *
 * S:            Symbolic analysis struct for QR factorization.
 *               On input it contains information that is not used in this
 *               function such as the row/column permutation
 *               On output it contains the number of nonzeros in R.
 */

SPEX_info spex_qr_pre_factor
(
    // Output
    SPEX_matrix *R_handle,        // On output: partial L matrix
                                  // On input: undefined
    // Input
    //int64_t *xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix A,          // Input Matrix
    const SPEX_symbolic_analysis S  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
)
{

    // All inputs have been checked by the caller, thus asserts are used here
    // as a reminder of the expected data types
    SPEX_info info;
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);

    //int64_t  top, k, j, jnew, n = A->n, p = 0;
    int64_t *w, *s, *leftmost;
    int64_t top, k, len, i, p, n = A->n, m2=n, m=A->m, rnz;
    //int64_t *c = NULL;
    SPEX_matrix R = NULL;
    ASSERT(n >= 0);

    //--------------------------------------------------------------------------
    // Declare memory for L and c
    //--------------------------------------------------------------------------

    // Allocate L
    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, S->unz,
        false, false, NULL));

    // Allocate c
    //c = (int64_t*) SPEX_malloc(n* sizeof (int64_t));

    w = (int64_t*) SPEX_malloc((n+m2)* sizeof (int64_t));
    leftmost = (int64_t*) SPEX_malloc(m* sizeof (int64_t));
    s = w + n ;
    if (!w)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // Set the column pointers of L and c
    for (k = 0; k < n; k++)
    {
        R->p[k] = S->cp[k+1];
    }

    R->i[0] = 0;
    //c[0]++;


    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */

    for (i = 0 ; i < m ; i++) leftmost [i] = -1 ;
    for (k = n-1 ; k >= 0 ; k--)
    {
        for (p = A->p [k] ; p < A->p [k+1] ; p++)
        {
            leftmost [A->i [p]] = k ;         /* leftmost[i] = min(find(A(i,:)))*/
        }
    }
    //--------------------------------------------------------------------------
    // Iterations 1:n-1
    //--------------------------------------------------------------------------
    rnz = 0 ;
    for (k = 0; k < n; k++)
    {
        w [k] = k ;  
        top = n ;

        for (p = A->p [k] ; p < A->p [k+1] ; p++)   /* find R(:,k) pattern */
        {
            i = leftmost [A->i [p]] ;         /* i = min(find(A(i,q))) */
            for (len = 0 ; w [i] != k ; i = S->parent [i]) /* traverse up to k */
            {
                s [len++] = i ;
                w [i] = k ;
            }
            while (len > 0) s [--top] = s [--len] ; /* push path on stack */
        }
        //from here to 
        i = A->i [p] ;             /* i = permuted row of A(:,col) */
        if (i > k && w [i] < k)         /* pattern of V(:,k) = x (k+1:m) */
        {
            w [i] = k ;
        }
        //here is prob not needed

        for (p = top ; p < n ; p++) /* for each i in pattern of R(:,k) */
        {
            i = s [p] ;                     /* R(i,k) is nonzero */
            R->i [rnz++] = i ;                  /* R(i,k) = x(i) */
        }
        R->i [rnz++] = k ;                     /* R(k,k) */
    }
    // Finalize R->p
    R->p[n] = S->unz = rnz;
    (*R_handle) = R;

    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
