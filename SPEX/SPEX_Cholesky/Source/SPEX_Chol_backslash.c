//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_backslash: solve Ax=b, returning solution as desired data\
//                                type
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* Purpose: This code utilizes the SPEX Cholesky factorization to exactly solve
 *          the linear system Ax = b.
 *
 * Input/Output arguments:
 *
 * x_handle:    A pointer to the solution of the linear system. The output is
 *              allowed to be returned in either double precision, mpfr_t, or
 *              rational mpq_t
 *
 * type:        Data structure of output desired. Must be either SPEX_MPQ,
 *              SPEX_FP64, or SPEX_MPFR
 *
 * A:           User's input matrix. It must be populated prior to calling this
 *              function.
 *
 * b:           Collection of right hand side vectors. Must be populated prior
 *              to factorization.
 *
 * option:      Struct containing various command parameters for the
 *              factorization. If NULL on input, default values are used.
 */

//TODO fix this everywhere else!! CONSISTENCY it should be SPEX_FREE_WORKSPACE
// functions that allocate for themselves is WORK
// functions that allocate "to return" is SPEX_FREE_ALLOCATION (instead of FREE_ALL) (to free workspce and allocations that were going to be returned)
// change to macro in util (and this will be propagated throughout)
# define SPEX_FREE_WORKSPACE       \
    SPEX_factorization_free(&F, option);     \
    //SPEX_symbolic_analysis_free (&S, option);

# define SPEX_FREE_ALLOCATION     \
    SPEX_FREE_WORKSPACE           \
    SPEX_matrix_free(&x, NULL);   \

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_backslash
(
    // Output
    SPEX_matrix** x_handle,       // On input: undefined. 
                                  // On output: final solution vector
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_MPQ, SPEX_MPFR, or SPEX_FP64
    const SPEX_matrix* A,         // Input matrix
    const SPEX_matrix* b,         // Right hand side vector(s)
    const SPEX_options* option    // Command options
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

    // x can't be NULL
    //TODO TIM: isn't the below condition always true? (same question as in the SPEX_analysis_free)
    //If it is always true, then it would be impossible to check-cover the guts of the if...
    if (!x_handle)
    {
        return SPEX_INCORRECT_INPUT;
    }
    /*if (!(*x_handle) || !A || !b || !option) //TOASK it gives me error when using this :/
    {
        return SPEX_INCORRECT_INPUT;
    }*/
    if (!A || !b || !option)
    {
        return SPEX_INCORRECT_INPUT;
    }

    // type must be acceptable
    if (type != SPEX_MPQ && type != SPEX_FP64 && type != SPEX_MPFR)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(b->type == SPEX_MPZ);
    ASSERT(b->kind == SPEX_DENSE);

    ASSERT(A->n >= 0);
    ASSERT(A->m >= 0);
    ASSERT(A->n == A->m);
    if (A->n == 0 || A->m == 0 || A->n != A->m)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    // Declare memory
    SPEX_symbolic_analysis *S = NULL;
    SPEX_factorization *F = NULL ;
    SPEX_matrix *x = NULL;
    SPEX_matrix* PAP = NULL;
    
    
    //--------------------------------------------------------------------------
    // Symbolic Analysis: obtain the row/column ordering of A
    //--------------------------------------------------------------------------
    SPEX_CHECK(SPEX_Chol_preorder(&S, A, option));
    
    //--------------------------------------------------------------------------
    // Determine if A is indeed symmetric. If so, we try Cholesky.
    // The symmetry check here checks both the nonzero pattern and values.
    // In addition, the symmetry check also checks that no diagonal entry is zero;
    // as otherwise this indicates that the matrix is not SPD (even if symmetric)
    // If the symmetry check fails, the appropriate error code is returned
    //--------------------------------------------------------------------------


    // TODO: Change this determine symmetry routine to check the diagonals as well
    SPEX_CHECK(SPEX_determine_symmetry( (SPEX_matrix*) A, 1, option));

    //--------------------------------------------------------------------------
    // Permute matrix A, that is apply the row/column ordering from the 
    // symbolic analysis step to get the permuted matrix PAP.
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_Chol_permute_A(&PAP, A, S));

    //--------------------------------------------------------------------------
    // SPEX Chol Factorization: Perform the REF Cholesky factorization of 
    // A. By default, up-looking Cholesky factorization is done; however,
    // the left looking factorization is done if option->algo=SPEX_CHOL_LEFT
    //-------------------------------------------------------------------------- 

    //TODO:The functions Factor and solve (below) should expect option to be const; DONE
    //TODO: After those changes are done, then the casting (SPEX_options*) needs to be removed from here. DONE
    SPEX_CHECK(SPEX_Chol_Factor(&F, S,PAP, option));

    //--------------------------------------------------------------------------
    // Solve: Solve Ax = b using the REF Cholesky factorization. That is,
    // apply the factorization LDL' = PAP' to solve the linear system LDL'x = b.
    // At the conclusion of the solve routines, x is the exact solution of 
    // Ax = b stored as a set of numerators and denominators (mpq_t)
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_Chol_Solve(&x, F, b, option));

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
        SPEX_CHECK(SPEX_matrix_copy(&x2, SPEX_DENSE, type, x, option));
        (*x_handle) = x2;
        SPEX_matrix_free (&x, NULL);
    }

    //--------------------------------------------------------------------------
    // Free all workspace and return success
    //--------------------------------------------------------------------------
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
