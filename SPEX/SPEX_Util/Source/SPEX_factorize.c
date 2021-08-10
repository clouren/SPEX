//------------------------------------------------------------------------------
// SPEX_Util/SPEX_factorize: exact sparse factorization based on the available
// symbolic analysis and its kind, which should be LU, Cholesky or QR.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs the SPEX factorization based on given
 * symbolic analysis. For detailed infomation that can be obtained from the
 * SPEX_factorization, please refer to the user manual for the corresponding
 * data structure.
 *
 * F: undefined on input, created on output
 *
 * A: input only, not modified
 * S: input only, not modified
 * option: input only, not modified
 */


#include "spex_util_internal.h"

SPEX_info SPEX_factorize
(
    SPEX_factorization **F_handle, // The resulted factorization as specified
    const SPEX_matrix *A, // Input matrix
    const SPEX_factorization *S, // symbolic analysis
    const SPEX_options *option  // Control parameters
)
{
    if (F_handle == NULL ||
        (S->kind != SPEX_LU_ANALYSIS && S->kind != SPEX_CHOLESKY_ANALYSIS &&
         S->kind != SPEX_QR_ANALYSIS))
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*F_handle) = NULL;
    SPEX_info info;

    // allocate memory space for SPEX_factorization struct
    (*F_handle) = (SPEX_factorization*) SPEX_calloc(1,
        sizeof(SPEX_factorization));
    if ((*F_handle) == NULL)
    {
        return SPEX_OUT_OF_MEMORY;
    }

    if (S->kind == SPEX_LU_ANALYSIS)
    {
        // perform spex left-looking lu factorization
        info = spex_left_lu_factorize(&((*F_handle)->L), &((*F_handle)->U),
            &((*F_handle)->rhos), &((*F_handle)->Pinv_perm),
            &((*F_handle)->P_perm), A, S->Q_perm, S->lnz, S->unz, option);
        if (info != SPEX_OK)
        {
            SPEX_factorization_free(F_handle);
            return info;
        }

        // copy column permutation from symbolic analysis to factorization
        (*F_handle)->Q_perm = (int64_t*) SPEX_malloc(A->n * sizeof(int64_t));
        if ((*F_handle)->Q_perm == NULL)
        {
            SPEX_factorization_free(F_handle);
            return SPEX_OUT_OF_MEMORY;
        }
        memcpy((*F_handle)->Q_perm, S->Q_perm, A->n * sizeof(int64_t));

        // set the kind for factorization
        (*F_handle)->kind = SPEX_LU_FACTORIZATION;
    }
    else if (S->kind == SPEX_CHOLESKY_ANALYSIS)
    {
        // perform spex Cholesky factorization
        info = SPEX_Chol_Factor(&((*F_handle)->L), &((*F_handle)->rhos),
            S /*TODO write a internal function to accept SPEX_factorization?*/,
            A, left /*TODO to be in the option*/, option);
        if (info != SPEX_OK)
        {
            SPEX_factorization_free(F_handle);
            return info;
        }

        // copy column permutation and its inverse from symbolic analysis
        // to factorization
        (*F_handle)->P_perm = (int64_t*) SPEX_malloc(A->n * sizeof(int64_t));
        (*F_handle)->Pinv_perm = (int64_t*) SPEX_malloc(A->n * sizeof(int64_t));
        if ((*F_handle)->P_perm == NULL || (*F_handle)->Pinv_perm == NULL)
        {
            SPEX_factorization_free(F_handle);
            return SPEX_OUT_OF_MEMORY;
        }
        memcpy((*F_handle)->P_perm, S->P_perm, A->n * sizeof(int64_t));
        memcpy((*F_handle)->Pinv_perm, S->Pinv_perm, A->n * sizeof(int64_t));

        // set the kind for factorization
        (*F_handle)->kind = SPEX_CHOLESKY_FACTORIZATION;
    }
    else // TODO QR
    {
        // set the kind for factorization
        (*F_handle)->kind = SPEX_QR_FACTORIZATION;
    }
    
    // copy the scale of A to factorization
    info = SPEX_mpq_init((*F_handle)->scale_for_A);
    if (info != SPEX_OK)
    {
        SPEX_factorization_free(F_handle);
        return info;
    }
    info = SPEX_mpq_set ((*F_handle)->scale_for_A, A->scale);
    if (info != SPEX_OK)
    {
        SPEX_factorization_free(F_handle);
        return info;
    }
    return SPEX_OK;
}
