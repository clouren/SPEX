//------------------------------------------------------------------------------
// SPEX_Util/SPEX_solve: exactly the given sparse linear system Ax=b using
// available factorization.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function exactly solve the given sparse linear system Ax=b
 * using the provided factorization, which could be LU, Cholesky or QR, and
 * either updatable or not updatable.
 *
 * x: undefined on input, created on output
 *
 * b: input only, not modified
 * F: input only, not modified
 * option: input only, not modified
 */


#include "spex_util_internal.h"

SPEX_info SPEX_solve
(
    SPEX_matrix **x_handle, // the solution to the system
    const SPEX_matrix *b,   // the right-hand-side vector
    const SPEX_factorization *F, // factorization of A
    const SPEX_options *option  // Control parameters
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    SPEX_info info;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    if (!x_handle || !F ||
        (F->kind != SPEX_LU_FACTORIZATION &&
         F->kind != SPEX_CHOLESKY_FACTORIZATION &&
         F->kind != SPEX_QR_FACTORIZATION))
    {
        return SPEX_INCORRECT_INPUT;
    }
    SPEX_REQUIRE (b,       SPEX_DENSE, SPEX_MPZ) ;
    SPEX_REQUIRE (F->rhos, SPEX_DENSE, SPEX_MPZ) ;
    int64_t n = b->m;

    (*x_handle) = NULL;

    //--------------------------------------------------------------------------
    // perform symbolic analysis for specified kind
    //--------------------------------------------------------------------------
    if (F->kind == SPEX_LU_SYMBOLIC_ANALYSIS)
    {
        SPEX_REQUIRE (F->L,    SPEX_CSC,   SPEX_MPZ) ;
        SPEX_REQUIRE (F->U,    SPEX_CSC,   SPEX_MPZ) ;
        if (!(F->Q_perm) || !(F->Pinv_perm) || F->L->m != n || F->L->n != n ||
            F->U->n != n || F->U->m != n || F->rhos->m != n)
        {
            return SPEX_INCORRECT_INPUT;
        }

        // TODO if (F->L->kind == SPEX_DYNAMIC_CSC)
        // 
        info = spex_left_lu_solve(x_handle, b, F, option);
        if (info != SPEX_OK)        {  return info;    }
    }
    else if (F->kind == SPEX_CHOLESKY_SYMBOLIC_ANALYSIS)
    {
        SPEX_REQUIRE (F->L,    SPEX_CSC,   SPEX_MPZ) ;
        if (!(F->P_perm) || !(F->Pinv_perm) || F->L->m != n || F->L->n != n ||
            F->rhos->m != n)
        {
            return SPEX_INCORRECT_INPUT;
        }

        // TODO: make this function internal?
        // TODO remove input checking that has been done here
        // TODO make sure the output solution is in mpz_t with pending scale
        // stored in x->scale, refer to spex_left_lu_solve
        // TODO modify the function interface as following, PAP is used for
        // PAP->scale which is F->scale_for_A; A is not needed since the
        // solution check will be done seperately; other components are in F.
        info = SPEX_Chol_Solve(x_handle, b, F, option)
    }
    else // TODO QR
    {
    }

    return SPEX_OK;
}
