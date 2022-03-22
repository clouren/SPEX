//------------------------------------------------------------------------------
// SPEX_Left_LU/SPEX_Left_LU_solve: exact solution of Ax=b
//------------------------------------------------------------------------------

// SPEX_Left_LU: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system LD^(-1)U x = b. This is an
 * internal function, please refer to SPEX_Left_LU_solve for user-callable. It
 * essnetially serves as a wrapper for all forward and backward substitution
 * routines. This function always returns the solution matrix x as a mpz_t
 * matrix. If a user desires to have double or mpfr output, they must create
 * a matrix copy.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Unitialized on input.
 *           on output, contains the exact rational solution of the system
 *
 * b:        Set of RHS vectors
 *
 * F:        LU factorization of A. Unmodified on input/output
 *
 * option:   command options
 */

#define SPEX_FREE_WORKSPACE             \
    SPEX_matrix_free (&b2, NULL) ;      \

#define SPEX_FREE_ALLOCATION            \
    SPEX_FREE_WORKSPACE                 \
    SPEX_matrix_free (&x, NULL) ;

#include "spex_left_lu_internal.h"

SPEX_info SPEX_Left_LU_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_matrix **x_handle,  // rational solution to the system
    // input:
    const SPEX_factorization* F, // LU factorization
    const SPEX_matrix *b,   // right hand side vector
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    SPEX_REQUIRE (b,    SPEX_DENSE, SPEX_MPZ) ;

    if (!x_handle || !F)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    SPEX_REQUIRE (F->L,    SPEX_CSC,   SPEX_MPZ) ;
    SPEX_REQUIRE (F->U,    SPEX_CSC,   SPEX_MPZ) ;
    SPEX_REQUIRE (F->rhos, SPEX_DENSE, SPEX_MPZ) ;
    
    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    *x_handle = NULL;
    int64_t n = F->L->n;

    SPEX_matrix *x = NULL;   // final solution
    SPEX_matrix *b2 = NULL;  // permuted b

    //--------------------------------------------------------------------------
    // b2 (Pinv_perm) = b
    //--------------------------------------------------------------------------

    SPEX_CHECK (spex_left_lu_permute_b (&b2, b, F->Pinv_perm, option)) ;

    //--------------------------------------------------------------------------
    // b2 = L\b2, via forward substitution
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_left_lu_forward_sub(F->L, b2, F->rhos));

    //--------------------------------------------------------------------------
    // b2 = b2 * det, where det=rhos[n-1]
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_matrix_mul(b2, F->rhos->x.mpz[n-1], option));

    //--------------------------------------------------------------------------
    // b2 = U\b2, via back substitution
    //--------------------------------------------------------------------------
    SPEX_CHECK(spex_left_lu_back_sub(F->U, b2));

    //--------------------------------------------------------------------------
    // x = Q*b2
    //--------------------------------------------------------------------------

    // TODO rename the following function?
    // it would be more efficient to just perform swapping
    SPEX_CHECK (spex_left_lu_permute_b (&x, b2, F->Q_perm, option)) ;
    SPEX_matrix_free (&b2, NULL) ;

    //--------------------------------------------------------------------------
    // Update the scale properly
    //--------------------------------------------------------------------------
    // set x->scale = b->scale * rhos[n-1] / A_scale
    SPEX_CHECK(SPEX_mpq_set(x->scale, b->scale));
    SPEX_CHECK(SPEX_mpz_mul(SPEX_MPQ_NUM(x->scale),
                            SPEX_MPQ_NUM(x->scale), F->rhos->x.mpz[n-1]));
    SPEX_CHECK(SPEX_mpq_canonicalize(x->scale));
    SPEX_CHECK(SPEX_mpq_div(x->scale, x->scale, F->scale_for_A));

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORKSPACE ;
    (*x_handle) = x ;
    return (SPEX_OK) ;
}

