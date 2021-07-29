//------------------------------------------------------------------------------
// SPEX_Backslash/SPEX_Backslash.c: Solve a system Ax=b
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* Purpose: Exactly solve sparse linear systems using SPEX Software package
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

#include "SPEX_Backslash.h"

SPEX_info SPEX_Backslash
(
    // Output
    SPEX_matrix **X_handle,       // Final solution vector
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_MPQ, SPEX_MPFR, or SPEX_FP64
    const SPEX_matrix *A,         // Input matrix
    const SPEX_matrix *b,         // Right hand side vector(s)
    const SPEX_options* option    // Command options
)
{
    SPEX_info info;
    // Check inputs
    if ( !spex_initialized () ) return (SPEX_PANIC);
    
    if (X_handle == NULL)
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*X_handle) = NULL;
    
    if (type != SPEX_MPQ && type != SPEX_FP64 && type != SPEX_MPFR)
    {
        return SPEX_INCORRECT_INPUT;
    }

    // A must be CSC and MPZ
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_CSC);

    // b must be dense and MPZ
    ASSERT( b->type == SPEX_MPZ);
    ASSERT( b->kind == SPEX_DENSE);
    
    // A must be a square matrix
    ASSERT(A->m == A->n);
    
    // Declare output
    SPEX_matrix* x = NULL;
    
    // Determine if A is symmetric by checking both the pattern
    // and values. The output of this function is either:
    // SPEX_OK:          Matrix is symmetric, try Cholesky
    // SPEX_UNSYMMETRIC: Matrix is unsymmetric, try LU
    // Other error code: Some error occured. Return the error
    info = SPEX_determine_symmetry( (SPEX_matrix*) A, 1);

    if (info == SPEX_OK)
    {
        // A was classified to be symmetric. Attempt a Cholesky
        // factorization of A
        
        // Since Cholesky is occuring, we update the option
        // struct to do AMD and diagonal pivoting
        spex_backslash_set_defaults(option, false);
        
        // Try SPEX Cholesky. The output for this function
        // is either:
        // SPEX_OK: Cholesky success, x is the exact solution
        // SPEX_SINGULAR: Cholesky factorization failed. This means
        //                either A is singular or simply that it's
        //                not SPD. In this case, we try LU
        // Other error code: Some error. Return the error code and exit
        info = SPEX_Chol_backslash(&x, type, A, b, option);
        if (info == SPEX_OK)
        {
            // Cholesky was successful. Set X_handle = x
            (*X_handle) = x;
        }
        else if (info == SPEX_SINGULAR)
        {
            // Cholesky factorization failed. Must try
            // LU factorization now
            
            // Since LU is occuring we update the options struct
            // to use COLAMD and tolerance based pivoting.
            spex_backslash_set_defaults(option, true);
            
            // The LU factorization can return either:
            // SPEX_OK: LU success, x is the exact solution
            // Other error code: Some error. Return the error
            //                   code and exit
            info = SPEX_Left_LU_backslash(&x, type, A, b, option);
            if (info == SPEX_OK)
            {
                // LU success, set X_handle = x
                (*X_handle) = x;
            }
            else
            {
                // Both Cholesky and LU have failed, info contains
                // the problem, most likely that A is singular
                return info;
            }
        }
        else
        {
            // Cholesky failed, but not due to a SPEX_SINGULAR 
            // error code. Most likely invalid input or out of 
            // memory condition.
            return info;
        }
    }
    else if (info == SPEX_UNSYMMETRIC)
    {
        // A is classifed as unsymmetric, so an LU factorization
        // will be used.
        
        
        // Since LU is occuring we update the options struct
        // to use COLAMD and tolerance based pivoting.
        spex_backslash_set_defaults(option, false);
        
        // The LU factorization can return either:
        // SPEX_OK: LU success, x is the exact solution
        // Other error code: Some error. Return the error
        //                   code and exit
        info = SPEX_Left_LU_backslash(&x, type, A, b, option);
        if (info == SPEX_OK)
        {
            // LU factorization was successful. Set X_handle = x
            (*X_handle) = x;
        }
        else
        {
            // LU factorization failed. Return the error code
            return info;
        }
    }
    else
    {
        // An error occured during the classification of A
        // Return the error code and stop
        return info;
    }
    
    // At this point, x is the exact solution of Ax = b and is
    // stored in the user desired type. Now, we exit.
    return SPEX_OK;
}
