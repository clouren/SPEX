//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_solve: find the exact solution for Ax=b with the
// the updatable LU factorizaiton of A.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system L(P,:)D^(-1)U(:,Q) x = b. It
 * essnetially serves as a wrapper for all forward and backward substitution
 * routines. This function always returns the solution matrix x as a mpz_t
 * matrix with additional pending scale factor. If a user desires to have
 * each value for each entries with pending scale factor applied, simply
 * compute x->v[j]->x[i]/x->scale and convert to double or mpfr output as
 * desired.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Memory space will be allocated
 *           for x_handle to store the exact solution of the system
 *
 * b:        Set of RHS vectors
 *
 * F:        SPEX LU factorization in the updatable format
 *
 * option:   command options
 */

#include "spex_update_internal.h"

SPEX_info SPEX_Update_solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    const SPEX_matrix *b,   // a m*n dense matrix contains the right-hand-side
                            // vector
    SPEX_factorization *F,// The SPEX LU factorization in updatable form
    const SPEX_options* option // Command options
)
{
    SPEX_info info ;
    SPEX_CHECK(spex_update_solve_internal(x_handle, b, F, false, option));
 }
