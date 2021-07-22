//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Factor: Integer preserving Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"  

/* Purpose: This function performs the integer preserving Cholesky factorization.
 * It allows either the left-looking or up-looking integer-preserving Cholesky factorization.
 * In order to compute the L matrix, it performs n iterations of a sparse REF symmetric
 * triangular solve function. The overall factorization is PAP' = LDL'
 * 
 * Input arguments:
 * 
 * L_handle:    A handle to the L matrix. Null on input. On output, contains a pointer to the 
 *              L matrix
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. On output,
 *              contains the elimination tree and estimate number of nonzeros in L.
 * 
 * rhos_handle: A handle to the sequence of pivots. NULL on input. On output, contains a pointer
 *              to the pivots matrix.
 *
 * A:           The user's permuted input matrix
 * 
 * left:        A boolean parameter which tells the function whether it is performing a left-looking
 *              or up-looking factorization. If this bool is true, a left-looking factorization
 *              is done, otherwise the up-looking factorization is done.
 * 
 * option:      Command options
 * 
 */
SPEX_info SPEX_Chol_Factor      
(
    // Output
    SPEX_matrix** L_handle,     // Lower triangular matrix. NULL on input
    SPEX_matrix** rhos_handle, // Sequence of pivots. NULL on input.
    SPEX_Chol_analysis* S,     // Symbolic analysis struct that contains elimination tree of A, column pointers of L, 
                                //exact number of nonzeros of L and permutation used
    // Input
    const SPEX_matrix* A,      // Matrix to be factored   
    bool left,                 // Set to true if performing a left-looking factorization; 
                               //otherwise perform an up-looking factorization.
    const SPEX_options* option //command options
)
{
    SPEX_info info;
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!L_handle || !S || !rhos_handle || !option )
    {
        return SPEX_INCORRECT_INPUT;
    }
    ASSERT(*L_handle==NULL);
    ASSERT(*rhos_handle==NULL);

    int64_t anz;
    SPEX_CHECK(SPEX_matrix_nnz (&anz, A, option)) ;
    
    if (anz < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*rhos_handle) = NULL ;
        
    //--------------------------------------------------------------------------
    // Declare outputs 
    //--------------------------------------------------------------------------
    
    SPEX_matrix *L = NULL ;
    SPEX_matrix *rhos = NULL ;
    
    //--------------------------------------------------------------------------
    // Call factorization
    //--------------------------------------------------------------------------
    
    // TODO Put this in option
    if (left)
    {
        SPEX_CHECK(spex_Chol_Left_Factor(&L, &rhos, S, A, option));
    }
    else
    {
        SPEX_CHECK(spex_Chol_Up_Factor(&L, &rhos, S, A, option));
    }
    
    //--------------------------------------------------------------------------
    // Set outputs, return ok
    //--------------------------------------------------------------------------
    
    (*L_handle) = L;
    (*rhos_handle) = rhos;
    return SPEX_OK;
}
