

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

    if (!F || !(F->P_perm) || !(F->Pinv_perm) || !(F->L) ||
        F->L->type != SPEX_MPZ ||
        (F->L->kind != SPEX_CSC && F->L->kind != SPEX_DYNAMIC_CSC) ||
        (F->kind != SPEX_LU_FACTORIZATION &&
         F->kind != SPEX_CHOLESKY_FACTORIZATION) ||
        (F->kind == SPEX_LU_FACTORIZATION && (!(F->Q_perm) || !(F->U) ||
         (F->updatable && !(F->Qinv_perm)) || F->U->type != SPEX_MPZ ||
         (F->U->kind != SPEX_CSC && F->U->kind != SPEX_DYNAMIC_CSC))))
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

