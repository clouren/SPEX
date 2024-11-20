//------------------------------------------------------------------------------
// SPEX_Update/SPEX_update_tsolve: find the exact solution for Ax=b with the
// the updatable LU factorizaiton of A.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2023, Christopher Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system A^T*x = b via LU or Cholesky
 * factorization of A. It essentially serves as a wrapper for all forward and
 * backward substitution routines. This function always returns the solution
 * matrix x as a mpq_t matrix. If a user desires to have other entry type,
 * simple call SPEX_matrix_copy.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Memory space will be allocated
 *           for x_handle to store the exact solution of the system
 *
 * F:        SPEX LU or Cholesky factorization of A
 *
 * b:        Set of RHS vectors
 *
 * option:   command options
 */

#include "spex_update_internal.h"

SPEX_info SPEX_update_tsolve // solves A^T*x = b
(
    // Output
    SPEX_matrix *x_handle,  // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    SPEX_factorization F,   // The SPEX LU or Cholesky factorization of A
    const SPEX_matrix b,    // a m*n dense matrix contains the right-hand-side
                            // vector
    const SPEX_options option // Command options
)
{
    return spex_update_solve_internal(x_handle, F, b, true, option);
}
