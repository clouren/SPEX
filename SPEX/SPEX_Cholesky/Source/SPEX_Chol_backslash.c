//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_backslash: solve Ax=b, returning solution as desired data
//                                type
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* Purpose: This code utilizes the SPEX Cholesky factorization to exactly solve
 *          the linear system Ax = b.
 *
 * Input/Output arguments:
 *
 * x_handle:    A pointer to the solution of the linear system. The output 
 *              can be returned in double precision, 
 *              mpfr_t (user-specified precision floating point), or
 *              mpq_t (rational)
 *
 * type:        Type of output desired. 
 *              Must be SPEX_MPQ, SPEX_FP64 or SPEX_MPFR
 *
 * A:           User's input matrix.
 *              Must be populated prior to calling this function.
 *
 * b:           Collection of right hand side vector(s).
 *              Must be populated prior to calling this function.
 *
 * option:      Struct containing various command parameters for the
 *              factorization. If NULL on input, default values are used.
 */

#define SPEX_FREE_WORKSPACE                 \
{                                            \
    SPEX_factorization_free(&F, option);     \
    SPEX_symbolic_analysis_free (&S, option);\
    SPEX_FREE (PAP->x.mpz);                  \
    SPEX_matrix_free(&PAP, NULL);            \
}  

#define SPEX_FREE_ALL            \
{                                 \
    SPEX_FREE_WORKSPACE           \
    SPEX_matrix_free(&x, NULL);   \
}

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_backslash
(
    // Output
    SPEX_matrix** x_handle,       // On input: undefined. 
                                  // On output: solution vector(s)
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_FP64, SPEX_MPFR, or SPEX_MPQ
    const SPEX_matrix* A,         // Input matrix. Must be SPEX_MPZ and SPEX_CSC
    const SPEX_matrix* b,         // Right hand side vector(s). Must be SPEX_MPZ and SPEX_DENSE
    const SPEX_options* option    // Command options (Default if NULL)
)
{

    SPEX_info info;
    // SPEX must be initialized
    if (!spex_initialized())
    {
        return SPEX_PANIC;
    }
    
    //-------------------------------------------------------------------------
    // check inputs
    //-------------------------------------------------------------------------

    // x, A, B can't be NULL
    if (!x_handle || !A || !b)
    {
        return SPEX_INCORRECT_INPUT;
    }

    // type must be acceptable
    if (type != SPEX_MPQ && type != SPEX_FP64 && type != SPEX_MPFR)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    // A must be the appropriate dimension
    if (A->n == 0 || A->m == 0 || A->n != A->m)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    // Make sure A and b are in the correct format
    if (A->type != SPEX_MPZ || A->kind != SPEX_CSC || b->type != SPEX_MPZ
        || b->kind != SPEX_DENSE)
    {
        return SPEX_INCORRECT_INPUT;
    }

    //TODO Make sure that all the functions below, handle a NULL oooptions
    // Is this done? If so please delete this
    
    // Declare memory
    SPEX_symbolic_analysis *S = NULL;
    SPEX_factorization *F = NULL ;
    SPEX_matrix *x = NULL;
    SPEX_matrix* PAP = NULL;
    
    
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
    // Factorization: Perform the REF Cholesky factorization of PAP.
    // By default, up-looking Cholesky factorization is done; however,
    // the left looking factorization is done if option->algo=SPEX_CHOL_LEFT
    //-------------------------------------------------------------------------- 

    SPEX_CHECK( spex_chol_factor(&F, S, PAP, option) );

    //--------------------------------------------------------------------------
    // Solve: Solve Ax = b using the REF Cholesky factorization. That is,
    // apply the factorization LDL' = PAP' to solve the linear system LDL'x = b.
    // At the conclusion of the solve function, x is the exact solution of 
    // Ax = b stored as a set of numerators and denominators (mpq_t)
    //--------------------------------------------------------------------------

    SPEX_CHECK( SPEX_Chol_solve(&x, F, b, option) );

    //--------------------------------------------------------------------------
    // At this point x is stored as mpq_t. If the user desires the output 
    // to be mpq_t we set x_handle = x. Otherwise, we create a copy, x2, of x
    // of the desired type. We then set x_handle = x2 and free x.
    //--------------------------------------------------------------------------

    if (type == SPEX_MPQ)
    {
        (*x_handle) = x;
    }
    else
    {
        SPEX_matrix* x2 = NULL;
        SPEX_CHECK( SPEX_matrix_copy(&x2, SPEX_DENSE, type, x, option) );
        (*x_handle) = x2;
        SPEX_matrix_free (&x, NULL);
    }


    //--------------------------------------------------------------------------
    // Free all workspace and return success
    //--------------------------------------------------------------------------
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
