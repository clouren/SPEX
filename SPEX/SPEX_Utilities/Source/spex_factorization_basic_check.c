//------------------------------------------------------------------------------
// SPEX_Utilities/spex_factorization_basic_check.c: perform basic check for a
// given factorization.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Chris Lourenco,
// Erick Moreno-Centeno.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_util_internal.h"

SPEX_info spex_factorization_basic_check
(
    SPEX_factorization F
)
{

    if (!spex_initialized()) {return SPEX_PANIC;}

    if (!F ||
        /* F can only be LU or Cholesky*/
        !(F->kind == SPEX_LU_FACTORIZATION ||
          F->kind == SPEX_CHOLESKY_FACTORIZATION) ||
        /* P_perm and Pinv_perm should exist for any kind*/ 
        !(F->P_perm) || !(F->Pinv_perm) ||
        /* L should exist and have MPZ entries*/
        !(F->L) || F->L->type != SPEX_MPZ ||
        /* L should be non-shallow csc for non-updatable or
                           dynamic csc for updatable*/
        !(
          (!(F->updatable) && F->L->kind == SPEX_CSC &&
                              !(F->L->p_shallow) && F->L->p &&
                              !(F->L->i_shallow) && F->L->i &&
                              !(F->L->x_shallow) && F->L->x.mpz     ) ||
          (  F->updatable  && F->L->kind == SPEX_DYNAMIC_CSC && F->L->v)
         ) ||
        /* if F is LU*/
        (F->kind == SPEX_LU_FACTORIZATION &&
        /* Q_perm should always exist and Qinv_perm should exist for updatable*/
         (!(F->Q_perm) || (F->updatable && !(F->Qinv_perm)) ||
        /*U should exist and have MPZ entries*/
          !(F->U) || F->U->type != SPEX_MPZ ||
        /* U should be non-shallow csc for non-updatable or
                           dynamic csc for updatable*/
          !(
            (!(F->updatable) && F->U->kind == SPEX_CSC &&
                                !(F->U->p_shallow) && F->U->p &&
                                !(F->U->i_shallow) && F->U->i &&
                                !(F->U->x_shallow) && F->U->x.mpz     ) ||
            (  F->updatable  && F->U->kind == SPEX_DYNAMIC_CSC && F->U->v)
           )                                                               )))
    {
        return SPEX_INCORRECT_INPUT;
    }

    return SPEX_OK ;
}

