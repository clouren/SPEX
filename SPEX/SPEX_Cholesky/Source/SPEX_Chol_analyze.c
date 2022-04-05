//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_analize: 
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: 
 *
 * Input/Output arguments:
 *
 */

# define SPEX_FREE_WORKSPACE                 \
{                                            \
    SPEX_matrix_free(&PAP, NULL); \
}

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_analyze
(
    // Output
    SPEX_symbolic_analysis** S_handle, // Symbolic analysis data structure 
    // Input
    const SPEX_matrix* A,         // Input matrix
    const SPEX_options* option    // Command options
)
{

    SPEX_info info;

    //TODO input check
    //SPEX_REQUIRE_kind()

    SPEX_matrix* PAP = NULL;
    SPEX_symbolic_analysis *S = NULL;


    //--------------------------------------------------------------------------
    // Preorder: obtain the row/column ordering of A
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_chol_preorder(&S, A, option));
    
    //--------------------------------------------------------------------------
    // Determine if A is indeed symmetric. If so, we try Cholesky.
    // The symmetry check here checks both the nonzero pattern and values.
    // In addition, the symmetry check also checks that no diagonal entry is zero;
    // as otherwise this indicates that the matrix is not SPD (even if symmetric)
    // If the symmetry check fails, the appropriate error code is returned
    // --------------------------------------------------------------------------

    SPEX_CHECK(SPEX_determine_symmetry( (SPEX_matrix*) A, 1, option));

    //--------------------------------------------------------------------------
    // Permute matrix A, that is apply the row/column ordering from the 
    // symbolic analysis step to get the permuted matrix PAP.
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_chol_permute_A(&PAP, A, false, S));

    //--------------------------------------------------------------------------
    // Symbolic Analysis: compute the elimination tree of A
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_chol_symbolic_analysis(S,PAP,option));

    //--------------------------------------------------------------------------
    // Free all workspace and return success
    //--------------------------------------------------------------------------
    (*S_handle) = S;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
