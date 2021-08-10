//------------------------------------------------------------------------------
// SPEX_Util/SPEX_analyze: symbolic ordering and analysis for specified
// factorization kind
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* SPEX_analyze performs the symbolic ordering and analysis for specified
 * factorization kind. The pre-ordering method is specified in the option.
 * For the pre-ordering for LU factorization, there are three options:
 * no ordering, COLAMD, and AMD.
 * TODO? For the pre-ordering for Cholesky factorization,
 *
 * Input/output arguments:
 *
 * S:       Symbolic analysis struct. Undefined on input; contains symbolic
 *          analysis results for specified factorization on output
 *
 * A:       Input matrix, unmodified on input/output
 *
 * kind:    The desired factorization kind
 *
 * option:  option->order tells the function which ordering scheme to use
 *
 */

#include "spex_util_internal.h"

SPEX_info SPEX_analyze
(
    SPEX_factorization **S_handle, // symbolic analysis
    const SPEX_matrix *A, // Input matrix
    const SPEX_factorization_kind kind, // LU, Cholesky or QR analysis
    const SPEX_options *option  // Control parameters
)
{
    // check inputs
    if (S_handle == NULL ||
        (kind != SPEX_LU_ANALYSIS && kind != SPEX_CHOLESKY_ANALYSIS &&
         kind != SPEX_QR_ANALYSIS))
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*S_handle) = NULL;
    SPEX_info info;

    // perform symbolic analysis for specified kind
    if (kind == SPEX_LU_ANALYSIS)
    {
        // perform symbolic analysis for LU factorization
        SPEX_LU_analysis *S = NULL;
        info = SPEX_LU_analyze(&S, A, option);
        if (info != SPEX_OK)
        {
            SPEX_analysis_free(&S, option);
            return info;
        }

        // allocate memory space for SPEX_factorization struct
        (*S_handle) = (SPEX_factorization*) SPEX_calloc(1,
            sizeof(SPEX_factorization));
        if ((*S_handle) == NULL)
        {
            SPEX_analysis_free(&S, option);
            return SPEX_OUT_OF_MEMORY;
        }

        // set SPEX_factorization kind
        (*S_handle)->kind   = kind;

        // copy symbolic analysis result to SPEX_factorization struct
        (*S_handle)->Q_perm = S->q;         S->q = NULL;

        (*S_handle)->lnz    = S->lnz;
        (*S_handle)->unz    = S->unz;

        // free SPEX_LU_analysis struct
        SPEX_LU_analysis_free(&S, option);
    }
    else if (kind == SPEX_CHOLESKY_ANALYSIS)
    {
        // perform symbolic analysis for Cholesky factorization
        SPEX_Chol_analysis *S = NULL;
        // TODO this is in SPEX_Cholesky
        info = SPEX_Chol_preorder(&S, A, option);
        if (info != SPEX_OK)
        {
            SPEX_Chol_analysis_free(&S);
            return info;
        }

        // allocate memory space for SPEX_factorization struct
        (*S_handle) = (SPEX_factorization*) SPEX_calloc(1,
            sizeof(SPEX_factorization));
        if ((*S_handle) == NULL)
        {
            SPEX_Chol_analysis_free(&S);
            return SPEX_OUT_OF_MEMORY;
        }

        // set SPEX_factorization kind
        (*S_handle)->kind      = kind;

        // copy symbolic analysis result to SPEX_factorization struct
        (*S_handle)->P_perm    = S->p;      S->p      = NULL;
        (*S_handle)->Pinv_perm = S->pinv;   S->pinv   = NULL;
        (*S_handle)->parent    = S->parent; S->parent = NULL;
        (*S_handle)->cp        = S->cp;     S->cp     = NULL;

        (*S_handle)->lnz       = S->lnz;

        // free SPEX_Chol_analysis struct
        SPEX_Chol_analysis_free(&S);
    }
    else // if (kind == SPEX_QR_ANALYSIS)
    {
        // TODO
    }

    return SPEX_OK;
}
