//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_analyze: Perform the symbolic analysis of A
//------------------------------------------------------------------------------

// TODO: Update this part for all makefiles (not from lab)

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: perform the symbolic analysis of A for the Cholesky factorization,
 * that is, preordering A, computing the elimination tree, getting the column
 * counts of A, setting the column pointers and exact number of non zeros of L.
 * 
 * Input arguments of the function:
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. 
 *              On input it's NULL
 *              On output it contains the row/column permutation, the elimination
 *              tree, and the number of nonzeros in L.
 * 
 * A:           User's input matrix (Must be SPEX_MPZ and SPEX_CSC)
 * 
 * option:      Command options (Default if NULL)
 * 
 */

# define SPEX_FREE_WORKSPACE      \
{                                 \
    SPEX_matrix_free(&PAP, NULL); \
}

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_analyze
(
    // Output
    SPEX_symbolic_analysis** S_handle, // Symbolic analysis data structure 
    // Input
    const SPEX_matrix* A,         // Input matrix. Must be SPEX_MPZ and SPEX_CSC
    const SPEX_options* option    // Command options (Default if NULL)
)
{

    SPEX_info info;
    // SPEX must be initialized
    if (!spex_initialized())
    {
        return SPEX_PANIC;
    }
    
    // Check inputs
    if ( !S_handle || !A)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    // SPEX must be CSC
    SPEX_REQUIRE_KIND(A, SPEX_CSC);

    // Declare permuted matrix and S
    SPEX_matrix* PAP = NULL;
    SPEX_symbolic_analysis *S = NULL;


    //--------------------------------------------------------------------------
    // Preorder: obtain the row/column ordering of A (Default is AMD)
    //--------------------------------------------------------------------------

    SPEX_CHECK( spex_chol_preorder(&S, A, option) );
    
    //--------------------------------------------------------------------------
    // Determine if A is indeed symmetric. If so, we try Cholesky.
    // This symmetry check checks for both the nonzero pattern and values.
    // In addition, the symmetry check also checks that no diagonal entry is zero;
    // as otherwise this indicates that the matrix is not SPD (even if symmetric)
    // If the symmetry check fails, the appropriate error code is returned
    //--------------------------------------------------------------------------

    SPEX_CHECK( SPEX_determine_symmetry((SPEX_matrix*)A, option) );

    //--------------------------------------------------------------------------
    // Permute matrix A, that is apply the row/column ordering from the 
    // symbolic analysis step to get the permuted matrix PAP.
    //--------------------------------------------------------------------------

    SPEX_CHECK( spex_chol_permute_A(&PAP, A, true, S) );

    //--------------------------------------------------------------------------
    // Symbolic Analysis: compute the elimination tree of PAP
    //--------------------------------------------------------------------------

    SPEX_CHECK( spex_chol_symbolic_analysis(S, PAP, option) );

    //--------------------------------------------------------------------------
    // Set output, free all workspace and return success
    //--------------------------------------------------------------------------

    (*S_handle) = S;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
