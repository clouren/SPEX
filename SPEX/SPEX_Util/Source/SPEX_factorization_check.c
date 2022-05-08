

#include "spex_util_internal.h"

SPEX_info SPEX_factorization_check
(
    SPEX_factorization *F, // The factorization to check
    const SPEX_options* option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!spex_initialized()) {return SPEX_PANIC;}

    //--------------------------------------------------------------------------
    // basic checks
    //--------------------------------------------------------------------------

    // TODO: maybe put this if test in a helper function
    // spex_factorization_basic_check (F)
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
                              !(F->L->p_shallow) && !(F->L->p) &&
                              !(F->L->i_shallow) && !(F->L->i) &&
                              !(F->L->x_shallow) && !(F->L->x)     ) ||
          (  F->updatable  && F->L->kind == SPEX_DYNAMIC_CSC && !(F->L->v))
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
                                !(F->U->p_shallow) && !(F->U->p) &&
                                !(F->U->i_shallow) && !(F->U->i) &&
                                !(F->U->x_shallow) && !(F->U->x)     ) ||
            (  F->updatable  && F->U->kind == SPEX_DYNAMIC_CSC && !(F->U->v))
           )                 )))
    {
        return SPEX_INCORRECT_INPUT;
    }


    //--------------------------------------------------------------------------
    // Check all F->L, etc matrices with SPEX_matrix_check
    //--------------------------------------------------------------------------

    SPEX_info info;

    info = SPEX_matrix_check (F->L, options) ;
    if (info != SPEX_OK) return (info) ;

    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        info = SPEX_matrix_check (F->U, options) ;
        if (info != SPEX_OK) return (info) ;
    }

    return (SPEX_OK) ;
}

