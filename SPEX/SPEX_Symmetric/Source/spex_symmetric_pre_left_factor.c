//------------------------------------------------------------------------------
// SPEX_Symmetric/spex_symmetric_pre_left_factor: symbolic left-looking
// Cholesky/LDL
// ------------------------------------------------------------------------------

// SPEX_Symmetric: (c) 2020-2025, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_FREE(cp);                  \
}

# define SPEX_FREE_ALL               \
{                                    \
    SPEX_FREE_WORKSPACE              \
    SPEX_matrix_free(&L, NULL);      \
}

#include "spex_symmetric_internal.h"

/* Purpose: This function performs a symbolic left-looking factorization.
 * On input, A is the matrix to be factored, parent contains the elimination
 * tree and S contains the row/column permutations and number of nonzeros in L.
 * On output, L_handle is allocated to contain the nonzero pattern of L and
 * memory for the values.
 *
 * Importantly, this function assumes that A has already been permuted.
 *
 * Input arguments of the function:
 *
 * L_handle:    A handle to the L matrix. Null on input.
 *              On output, contains a pointer to the partial L matrix.
 *
 * xi:          Workspace nonzero pattern vector. It stores the pattern of
 *              nonzeros of the kth column of L for the triangular solve.
 *
 * A:           The user's permuted input matrix
 *
 * S:            Symbolic analysis struct for Cholesky or LDL factorization.
 *               On input it contains information that is not used in this
 *               function such as the row/column permutation
 *               On output it contains the number of nonzeros in L.
 */

SPEX_info spex_symmetric_pre_left_factor
(
    // Output
    SPEX_matrix *L_handle,        // On output: partial L matrix
                                  // On input: undefined
    // Input
    int64_t *xi,                  // Workspace nonzero pattern vector
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

    int64_t top, n = A->n ;
    int64_t *cp = NULL;
    SPEX_matrix L = NULL;
    ASSERT(n >= 0);

    //--------------------------------------------------------------------------
    // Declare memory for L and cp
    //--------------------------------------------------------------------------

    // Allocate L
    SPEX_CHECK(SPEX_matrix_allocate(&L, SPEX_CSC, SPEX_MPZ, n, n, S->lnz,
        false, false, NULL));

    // Allocate cp
    cp = (int64_t*) SPEX_malloc(n* sizeof (int64_t));
    if (!cp)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // Set the column pointers of L and cp
    for (int64_t j = 0; j < n; j++)
    {
        L->p[j] = cp[j] = S->cp[j];
    }

    // cp[j] points to the first 'empty' entry in the jth column of L, where
    // the next entry L(k,j) must be placed.

    //--------------------------------------------------------------------------
    // Construct the pattern of L one row at a time
    //--------------------------------------------------------------------------

    // L(0,:): first row of L contains just the diagonal entry;
    // add this entry to L(:,0), the first column of L
    L->i[0] = 0;
    cp[0]++;

    // now handle rows 1 to n-1 of L
    for (int64_t k = 1; k < n; k++)
    {

        //----------------------------------------------------------------------
        // Obtain nonzero pattern of the kth row of L (L(k,0:k-1)) in xi[top..n]
        //----------------------------------------------------------------------

        SPEX_CHECK(spex_symmetric_ereach(&top, xi, A, k, S->parent, cp));

        //----------------------------------------------------------------------
        // Copy the entries in L(k,0:k-1) into their corresponding columns
        //----------------------------------------------------------------------

        for (int64_t px = top; px < n; px++)
        {
            // get L(k,j) from the stack xi
            int64_t j = xi[px];
            ASSERT (j < k) ;
            // place L(k,j) in the kth column of L
            int64_t p = cp[j]++;
            L->i[p] = k;
        }

        //----------------------------------------------------------------------
        // place L(k,k) in the kth column of L
        //----------------------------------------------------------------------

        int64_t p = cp[k]++;
        L->i[p] = k;
    }

    //--------------------------------------------------------------------------
    // Finalize L->p and return result
    //--------------------------------------------------------------------------

    L->p[n] = S->lnz;
    (*L_handle) = L;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
