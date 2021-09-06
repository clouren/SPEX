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
    // check inputs
    if (x_handle == NULL ||
        (F->kind != SPEX_LU_FACTORIZATION &&
         F->kind != SPEX_CHOLESKY_FACTORIZATION &&
         F->kind != SPEX_QR_FACTORIZATION))
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*x_handle) = NULL;
    SPEX_info info;

    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        // TODO if (F->L->kind == SPEX_DYNAMIC_CSC)
        // 
        info = spex_left_lu_solve(x_handle, b, F->scale_for_A, F->L, F->U,
            F->rhos, F->Q_perm, F->Pinv_perm, option);
        if (info != SPEX_OK)
        {
            SPEX_matrix_free(x_handle);
            return info;
        }
    }
    else if (F->kind == SPEX_CHOLESKY_FACTORIZATION)
    {
    }
    else // TODO QR
    {
    }

    return SPEX_OK;
}
