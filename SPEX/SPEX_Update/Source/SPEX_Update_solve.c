//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_solve: find the exact solution for Ax=b with the
// the updatable LU factorizaiton of A.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system A*x = b via the LU or
 * Cholesky factorization of A. It essnetially serves as a wrapper for all
 * forward and backward substitution routines. This function always returns the
 * solution matrix x as a mpq_t matrix. If a user desires to have other entry
 * type, simple call SPEX_matrix_copy.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Memory space will be allocated
 *           for x_handle to store the exact solution of the system
 *
 * b:        Set of RHS vectors
 *
 * F:        SPEX LU or Cholesky factorization
 *
 * option:   command options
 */

#include "spex_update_internal.h"

SPEX_info SPEX_Update_solve // solves Ax = b via LU or Cholesky factorization
(
    // Output
    SPEX_matrix **x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    SPEX_factorization *F,  // The SPEX LU or Cholesky factorization
    const SPEX_matrix *b,   // a m*n dense matrix contains the right-hand-side
                            // vector
    const SPEX_options* option // Command options
)
{
    return spex_update_solve_internal(x_handle, b, F, false, option);
}
