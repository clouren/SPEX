//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_factorize: Perform the SPEX Chol factorization of A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs the integer preserving Cholesky factorization.
 * First it permutes the input matrix according to the symbolic analysis.
 * Then it performs the left-looking or up-looking integer-preserving Cholesky
 * factorization. In order to compute the L matrix, it performs n iterations of
 * a sparse REF symmetric triangular solve function. The overall factorization
 * is PAP' = LDL'
 * 
 * 
 * Input arguments of the function:
 * 
 * F_handle:    Handle to the factorization struct. Null on input.
 *              On output, contains a pointer to the factorization.
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. 
 *              On input it contains the elimination tree and 
 *              the number of nonzeros in L.
 *
 * A:           User's input matrix. Must be SPEX_MPZ and SPEX_CSC.
 * 
 * option:      Command options. Default if NULL. Notably, option->chol_type indicates whether
 *              it is performing the default up-looking factorization (SPEX_CHOL_UP)
 *              or the left-looking factorization (SPEX_CHOL_LEFT).
 * 
 */

# define SPEX_FREE_WORKSPACE      \
{                                 \
    SPEX_FREE (PAP->x.mpz);       \
    SPEX_matrix_free(&PAP, NULL); \
}

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_factorize
(
    // Output
    SPEX_factorization **F_handle,  // Cholesky factorization struct
    //Input
    const SPEX_matrix* A,      // Matrix to be factored. Must be SPEX_MPZ and SPEX_CSC   
    const SPEX_symbolic_analysis* S,// Symbolic analysis struct containing the
                               // elimination tree of A, the column pointers of
                               // L, and the exact number of nonzeros of L.
    const SPEX_options* option //command options
                               // Notably, option->chol_type indicates whether
                               // CHOL_UP (default) or CHOL_LEFT is used.
)
{

    SPEX_info info;

    SPEX_matrix* PAP = NULL;
    SPEX_factorization *F = NULL ;

    //TODO check inputs with if
    ASSERT(S->kind==SPEX_CHOLESKY_FACTORIZATION);
    
    // TODO Copy comments from backslash
    //TODO copy space between comment to code from here

    //--------------------------------------------------------------------------
    // Numerically permute matrix A, that is apply the row/column ordering from 
    // the symbolic analysis step to get the permuted matrix PAP.
    //--------------------------------------------------------------------------
    SPEX_CHECK(spex_chol_permute_A(&PAP, A, true, S));

    //--------------------------------------------------------------------------
    // Factorization: Perform the REF Cholesky factorization of 
    // A. By default, up-looking Cholesky factorization is done; however,
    // the left looking factorization is done if option->algo=SPEX_CHOL_LEFT
    //-------------------------------------------------------------------------- 
    SPEX_CHECK(spex_chol_factor(&F, S,PAP, option));

    //--------------------------------------------------------------------------
    // Set F_handle = F, free all workspace and return success
    //--------------------------------------------------------------------------
    (*F_handle) = F ;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
