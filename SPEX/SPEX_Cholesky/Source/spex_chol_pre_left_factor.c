//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Pre_Left_Factor: Symbolic left-looking Chol factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


#include "spex_chol_internal.h"

/* Purpose: This function performs a symbolic left-looking factorization.
 * It allocates the memory for the L matrix and determines the full nonzero
 * pattern of L
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
 * S:            Symbolic analysis struct for Cholesky factorization. 
 *               On input it contains information that is not used in this 
 *               function such as the row/column permutation
 *               On output it contains the number of nonzeros in L.
 * 
 * c:            Column pointers
 * 
 */
SPEX_info spex_chol_pre_left_factor
(
    // Output
    SPEX_matrix** L_handle,       // On output: partial L matrix 
                                  // On input: undefined
    // Input
    int64_t* xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix* A,         // Input Matrix
    const SPEX_symbolic_analysis* S,  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
    int64_t* c                    // Column pointers
)
{
    // Input check
    SPEX_info info;
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);
    
    int64_t  top, k, j, jnew, n = A->n;
    ASSERT(n >= 0);
    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
       
    // Allocate L  
    SPEX_matrix* L = NULL;
    SPEX_matrix_allocate(&L, SPEX_CSC, SPEX_MPZ, n, n, S->lnz, false, false, NULL);
    for (k = 0; k < n; k++)
    {
        L->p[k] = (S->c)[k] = S->cp[k];
    }
        
    L->i[0] = 0;
    (S->c)[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        // Obtain nonzero pattern in xi[top..n]
        SPEX_CHECK(spex_chol_ereach(&top, xi, A, k, S->parent, (S->c)));
     
        //----------------------------------------------------------------------
        // Iterate accross the nonzeros in x
        //----------------------------------------------------------------------
        int64_t p = 0;
        for (j = top; j < n; j++)
        {
            jnew = xi[j];
            if (jnew == k) continue;
            p = (S->c)[jnew]++;
            // Place the i location of the L->nz nonzero
            L->i[p] = k;
        }
        p = (S->c)[k]++;
        L->i[p] = k;
    }
    // Finalize L->p
    L->p[n] = S->lnz;
    (*L_handle) = L;
    return SPEX_OK;
}
