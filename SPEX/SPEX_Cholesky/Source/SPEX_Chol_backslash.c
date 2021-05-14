//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_backslash: solve Ax=b, returning solution as desired data type
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------


/* Purpose: This code utilizes the SPEX Cholesky factorization to exactly solve the
 *          linear system Ax = b. This is essentially an exact version of
 *          MATLAB sparse backslash
 *
 * Input/Output arguments:
 *
 * X_handle:    A pointer to the solution of the linear system. The output is
 *              allowed to be returned in either double precision, mpfr_t, or
 *              rational mpq_t
 *
 * type:        Data structure of output desired. Must be either SPEX_MPQ,
 *              SPEX_FP64, or SPEX_MPFR
 *
 * A:           User's input matrix. It must be populated prior to calling this
 *              function.
 *
 * b:           Collection of right hand side vectors. Must be populated prior to
 *              factorization.
 *
 * option:      Struct containing various command parameters for the factorization. If
 *              NULL on input, default values are used.
 */

# define SPEX_FREE_WORK                 \
    SPEX_matrix_free(&L, NULL);         \
    SPEX_matrix_free(&A2,NULL);         \
    SPEX_matrix_free(&rhos, NULL);      \
    SPEX_Chol_analysis_free (&S);   \

# define SPEX_FREE_ALL              \
    SPEX_FREE_WORK                  \
    SPEX_matrix_free(&x, NULL);     \

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_backslash
(
    // Output
    SPEX_matrix** X_handle,       // Final solution vector
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_MPQ, SPEX_MPFR, or SPEX_FP64
    const SPEX_matrix* A,         // Input matrix
    const SPEX_matrix* b,         // Right hand side vector(s)
    const SPEX_options* option    // Command options
)
{
    //-------------------------------------------------------------------------
    // check inputs
    //-------------------------------------------------------------------------
    printf("chcek input\n");

    SPEX_info ok ;
    // SPEX must be initialized
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    ASSERT(*X_handle==NULL);

    // X can't be NULL
    if (X_handle == NULL)
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*X_handle) = NULL;

    // type must be acceptable
    if (type != SPEX_MPQ && type != SPEX_FP64 && type != SPEX_MPFR)
    {
        return SPEX_INCORRECT_INPUT;
    }

    SPEX_REQUIRE (A, SPEX_CSC,   SPEX_MPZ) ;
    SPEX_REQUIRE (b, SPEX_DENSE, SPEX_MPZ) ;

    SPEX_matrix *L = NULL ;
    SPEX_matrix *x = NULL;
    SPEX_matrix* A2 = NULL;
    SPEX_matrix *rhos = NULL ;
    SPEX_Chol_analysis *S = NULL;
    int64_t k, n = A->n, index;

    // n must be at least 0
    ASSERT(n >= 0);
    
    //--------------------------------------------------------------------------
    // Symbolic Analysis
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_Chol_preorder(&S, (SPEX_matrix*) A, option));
    
    //--------------------------------------------------------------------------
    // Determine if A is indeed symmetric. If so, we try Cholesky
    // --------------------------------------------------------------------------
    
    SPEX_CHECK(SPEX_determine_symmetry( (SPEX_matrix*) A, 1));    // Determine symmetry with nonzero pattern and values

    //--------------------------------------------------------------------------
    // Permute matrix A, that is set A2 = PAP'
    //--------------------------------------------------------------------------

    SPEX_CHECK( SPEX_Chol_permute_A(&A2, (SPEX_matrix*) A, S));

    //--------------------------------------------------------------------------
    // SPEX Chol Factorization
    //-------------------------------------------------------------------------- 
    
    SPEX_CHECK(SPEX_Chol_Factor(&L, &rhos, S,A2, 
                                 false,     // True = left, false = up
                                 (SPEX_options*) option));

    //--------------------------------------------------------------------------
    // Solve
    //--------------------------------------------------------------------------

    SPEX_CHECK (SPEX_Chol_Solve (&x, A2, (SPEX_matrix*) A, (SPEX_matrix*) b,
                                    rhos,
                                    L,
                                    S,
                                    (SPEX_options*) option));

    //--------------------------------------------------------------------------
    // Now, x contains the exact solution of the linear system in mpq_t
    // precision set the output.
    //--------------------------------------------------------------------------

    if (type == SPEX_MPQ)
    {
        (*X_handle) = x ;
    }
    else
    {
        SPEX_matrix* x2 = NULL ;
        SPEX_CHECK (SPEX_matrix_copy (&x2, SPEX_DENSE, type, x, option)) ;
        (*X_handle) = x2 ;
        SPEX_matrix_free (&x, NULL) ;
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    return (SPEX_OK) ;
}
