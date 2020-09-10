//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Pre_Left_Factor: Symbolic left-looking Chol
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------


#include "spex_chol_internal.h"

/* Purpose: This function performs a symbolic left-looking factorization
 * It allocates the memory for the L matrix and allocates the individual
 * entries in the matrix.
 */
SPEX_info spex_Chol_Pre_Left_Factor // pre-allocates the left-looking Chol factor
(
    SPEX_matrix* A,                 // Input matrix
    SPEX_matrix** L_handle,         // partial L matrix
    int64_t*  xi,                   // nonzero pattern vector
    int64_t*  parent,               // Elimination tree
    SPEX_Chol_analysis * S,         // stores nnz and elimination tree
    int64_t * c                     // Column pointers
)
{
    // Input check/
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
        top = spex_Chol_ereach(A, k, parent, xi, c);  // Obtain nonzero pattern in xi[top..n]
     
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
