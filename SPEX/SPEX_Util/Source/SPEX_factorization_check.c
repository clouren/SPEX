//------------------------------------------------------------------------------
// SPEX_Util/SPEX_factorization_check.c: check the correctness of a given
// factorization.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Chris Lourenco
// (US Naval Academy), Erick Moreno-Centeno, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_factorization_check checks all the followings for a give factorization:
// 1. if the required components exist;
//    (TODO use spex_factorization_basic_check?)
// 2. if sizes of different matrices match, i.e., L and U should be n*n, and
//    rhos should be n*1;
// 3. if L (and U if exists) is correct (using SPEX_matrix_check);
// 4. if L, U, and rhos have same pivot values, and when F is updatable, if L
//    (and U if exists) is of SPEX_DYNAMIC_CSC, and if the i-th pivot is the
//    first entry in the nonzero list of i-th vector of L (and U if exists);
// 5. if each permutation is reasonable, i.e., no duplicate, and in range of
//    [0,n), and if P_perm and Pinv_perm are mutually inverse vectors,
//    same applied to (Q_perm, Qinv_perm) if exists.


#include "spex_util_internal.h"

SPEX_info SPEX_factorization_check
(
    SPEX_factorization *F, // The factorization to check
    const SPEX_options* option
)
{

    if (!spex_initialized()) {return SPEX_PANIC;}

    //--------------------------------------------------------------------------
    // 1. check if the required components exist;
    //--------------------------------------------------------------------------

    SPEX_info info;

    // TODO: maybe put this if test in a helper function
    // spex_factorization_basic_check (F)
    if (!F ||
        /* F can only be LU or Cholesky*/
        !(F->kind == SPEX_LU_FACTORIZATION ||
          F->kind == SPEX_CHOLESKY_FACTORIZATION) ||
        /* rhos, P_perm and Pinv_perm should exist for any kind*/ 
        !(F->rhos) || !(F->P_perm) || !(F->Pinv_perm) ||
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
           )                 )))
    {
        return SPEX_INCORRECT_INPUT;
    }


    //--------------------------------------------------------------------------
    // 2. check if sizes of different matrices match, i.e., L and U should be
    //    n*n, and rhos should be n*1;
    //--------------------------------------------------------------------------

    int64_t i, j, lp, up, n = F->L->n;
    if (F->L->m != n || // L and U must be n-by-n
        (F->kind == SPEX_LU_FACTORIZATION && (F->U->m != n || F->U->n != n)) ||
        F->rhos->m != n || F->rhos->n != 1) // rhos must be n-by-1
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // 3. check if L, rhos (and U if exists) is correct
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_matrix_check (F->L, option)) ;
    SPEX_CHECK(SPEX_matrix_check (F->rhos, option)) ;

    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        SPEX_CHECK(SPEX_matrix_check (F->U, option)) ;
    }

    //--------------------------------------------------------------------------
    // 4. check if L, U, and rhos have same pivot values, and when F is
    //    updatable, if L (and U if exists) is of SPEX_DYNAMIC_CSC, and if the
    //    i-th pivot is the first entry in the nonzero list of i-th vector of L
    //    (and U if exists);
    //--------------------------------------------------------------------------
#if 0
    if (F->updatable)
    {
        for (int
    }
    else // non-updatable factorization
    {
        for (i = 0; i < n; i++)
        {
            for (lp = L->p[i]; lp < L->p[i+1]; lp++)
            {
                if (L->i[p]
            }
        }
    }
#endif

    //--------------------------------------------------------------------------
    // 5. check if each permutation is reasonable, i.e., no duplicate, and in
    //    range of [0,n), and if P_perm and Pinv_perm are mutually inverse
    //    vectors, same applied to (Q_perm, Qinv_perm) if exists.
    //--------------------------------------------------------------------------
    return (SPEX_OK) ;
}

