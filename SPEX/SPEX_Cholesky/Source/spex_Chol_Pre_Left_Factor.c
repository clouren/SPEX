//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Pre_Left_Factor: Symbolic left-looking Chol
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


#include "spex_chol_internal.h"

/* Purpose: This function performs a symbolic left-looking factorization
 * It allocates the memory for the L matrix and allocates the individual
 * entries in the matrix.
 */
SPEX_info spex_Chol_Pre_Left_Factor
(
    // Output
    SPEX_matrix** L_handle,       // On output: partial L matrix 
                                  // On input: undefined
    // Input
    int64_t* xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix* A,         // Input Matrix
    const int64_t* parent,        // Elimination tree
    const SPEX_Chol_analysis* S,  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, and the
                                  // row/coluimn permutation and its inverse  
    int64_t* c                    // Column pointers
)
{
    // TODO Remove this, make left-factor combine these two. (Chris)
    // Input check/
    SPEX_info info;
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!L_handle || !xi || !parent || !S || !c)
        return SPEX_INCORRECT_INPUT;
    
    int64_t  top, k, j, jnew, n = A->n;
    ASSERT(n >= 0);
    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
       
    // Allocate L  
    SPEX_matrix* L = NULL;
    SPEX_matrix_allocate(&L, SPEX_CSC, SPEX_MPZ, n, n, S->lnz, false, false, NULL);
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
        
    L->i[0] = 0;
    c[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        SPEX_CHECK(spex_Chol_ereach(&top, xi, A, k, parent, c));  // Obtain nonzero pattern in xi[top..n]
     
        //----------------------------------------------------------------------
        // Iterate accross the nonzeros in x
        //----------------------------------------------------------------------
        int64_t p = 0;
        for (j = top; j < n; j++)
        {
            jnew = xi[j];
            if (jnew == k) continue;
            p = c[jnew]++;
            // Place the i location of the L->nz nonzero
            L->i[p] = k;
        }
        p = c[k]++;
        L->i[p] = k;
    }
    // Finalize L->p
    L->p[n] = S->lnz;
    (*L_handle) = L;
    return SPEX_OK;
}
