//------------------------------------------------------------------------------
// SPEX_Backslash/SPEX_solve.c: Given a factorization, solve a system Ax=b
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2020-2025, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: Exactly solve a system of equations using a SPEX factorization.
 *          This function can be called after any of the spex factorizations
 *          and essentially serves as a wrapper for the other individual solves
 *          located in each folder.
 *
 * Input/Output arguments:
 *
 * x_handle:    Pointer to the exact solution of the system Ax = b
 *
 * F:           The factorization of A. Can be Cholesky, LDL, LU, or in the future QR
 *
 * b:           Collection of right hand side vectors.
 *
 * option:      Struct containing various command parameters for the
 *              factorization. If NULL on input, default values are used.
 */

#include "spex_util_internal.h"
#include "SPEX.h"

SPEX_info SPEX_solve
(
    // Output
    SPEX_matrix *x_handle,      // On input: undefined.
                                // On output: Rational solution (SPEX_MPQ)
                                // to the system.
    // input/output:
    SPEX_factorization F,       // The Cholesky or LDL factorization of A
    // input:
    const SPEX_matrix b,        // Right hand side vector
    const SPEX_options option   // command options
)
{

    SPEX_info info;
    // Check inputs
    if (!spex_initialized()) return SPEX_PANIC;

    (*x_handle) = NULL ;

    if (!F || !b) return SPEX_INCORRECT_INPUT;

    // b must be a dense mpz matrix
    SPEX_REQUIRE(b, SPEX_DENSE, SPEX_MPZ);

    // Ensure that F is properly structured
    SPEX_CHECK( SPEX_factorization_check(F, option));

    // The final solution
    SPEX_matrix x = NULL;

    // Determine the type of factorization used
    // This can be either LU, Cholesky, LDL, or QR
    SPEX_factorization_kind f_kind = F->kind;

    switch (f_kind)
    {
        // LU factorization was used. In this case, we call
        // LU forward and back solve to solve the system
        case SPEX_LU_FACTORIZATION:
            info = SPEX_lu_solve(&x, F, b, option);
            break ;

        // Cholesky factorization was used. In this case, we utilize
        // the Cholesky forward and back solve
        case SPEX_CHOLESKY_FACTORIZATION:
            info = SPEX_cholesky_solve(&x, F, b, option);
            break ;

        // LDL Factorization. Here we use the LDL forward and back solve
        case SPEX_LDL_FACTORIZATION:
            info = SPEX_ldl_solve(&x, F, b, option);
            break ;

        // QR factorization is requested. Not currently supported, thus
        // we return an error code
        case SPEX_QR_FACTORIZATION:
            info = SPEX_INCORRECT_INPUT;
            break ;
    }

    // x contains either the exact solution of the system or is NULL
    (*x_handle) = x;
    // returns SPEX_OK if the algorithm is successful or the appropriate error.
    return info;
}
