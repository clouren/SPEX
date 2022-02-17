//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_factorize: 
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

SPEX_info SPEX_Chol_factorize
(
    // Output
    SPEX_factorization **F_handle, // Cholesky factorization
    //Input
    const SPEX_symbolic_analysis* S,// Symbolic analysis struct containing the
                               // elimination tree of A, the column pointers of
                               // L, and the exact number of nonzeros of L.
    const SPEX_matrix* A,      // Matrix to be factored   
    const SPEX_options* option //command options
                               // Notably, option->chol_type indicates whether
                               // CHOL_UP (default) or CHOL_LEFT is used.
)
{

    SPEX_info info;

    SPEX_matrix* PAP = NULL;
    SPEX_factorization *F = NULL ;

    //--------------------------------------------------------------------------
    // Numerically permute matrix A, that is apply the row/column ordering from 
    // the symbolic analysis step to get the permuted matrix PAP.
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_chol_permute_A(&PAP, A, S));

    //--------------------------------------------------------------------------
    // Factorization: Perform the REF Cholesky factorization of 
    // A. By default, up-looking Cholesky factorization is done; however,
    // the left looking factorization is done if option->algo=SPEX_CHOL_LEFT
    //-------------------------------------------------------------------------- 

    SPEX_CHECK(spex_chol_factor(&F, S,PAP, option));

    //--------------------------------------------------------------------------
    // Free all workspace and return success
    //--------------------------------------------------------------------------
    (*F_handle) = F ;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}