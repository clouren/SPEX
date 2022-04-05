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
    const SPEX_symbolic_analysis *S, // symbolic analysis
    const SPEX_options *option  // Control parameters
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    SPEX_info info;
    if (!spex_initialized()) return SPEX_PANIC;
 
    SPEX_REQUIRE (A, SPEX_CSC, SPEX_MPZ) ;

    if (F_handle == NULL || S == NULL
        (S->kind != SPEX_LU_SYMBOLIC_ANALYSIS &&
         S->kind != SPEX_CHOLESKY_SYMBOLIC_ANALYSIS &&
         S->kind != SPEX_QR_SYMBOLIC_ANALYSIS))
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*F_handle) = NULL;

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    if (S->kind == SPEX_LU_SYMBOLIC_ANALYSIS)
    {
        // perform spex left-looking lu factorization
        info = spex_left_lu_factorize(F_handle, A, S, option);
        if (info != SPEX_OK)        {  return info;    }
    }
    else if (S->kind == SPEX_CHOLESKY_SYMBOLIC_ANALYSIS)
    {
        // perform spex Cholesky factorization
        // TODO: set the kind for factorization in the following function
        //(*F_handle)->kind = SPEX_CHOLESKY_FACTORIZATION;
        // TODO set F->scale_for_A
        // TODO remove input checking that has been done here
        // TODO this should be internal function?
        info = SPEX_Chol_Factor(F_handle,
            S /*TODO update the function to accept SPEX_symbolic_analysis?*/,
            A, left /*TODO to be in the option*/, option);
        if (info != SPEX_OK)        {  return info;    }

        // TODO this should be done inside SPEX_Chol_Factor
        // copy column permutation and its inverse from symbolic analysis
        // to factorization
        /*
        (*F_handle)->P_perm = (int64_t*) SPEX_malloc(A->n * sizeof(int64_t));
        (*F_handle)->Pinv_perm = (int64_t*) SPEX_malloc(A->n * sizeof(int64_t));
        if ((*F_handle)->P_perm == NULL || (*F_handle)->Pinv_perm == NULL)
        {
            SPEX_factorization_free(F_handle);
            return SPEX_OUT_OF_MEMORY;
        }
        memcpy((*F_handle)->P_perm, S->P_perm, A->n * sizeof(int64_t));
        memcpy((*F_handle)->Pinv_perm, S->Pinv_perm, A->n * sizeof(int64_t));
        */

        // TODO: should we add "parent" as a component of SPEX_factorization?
        // TODO: if so, copy this inside the above function
        // copy elimination tree of target matrix to factorization
        /*
        (*F_handle)->parent = (int64_t*) SPEX_malloc(A->n * sizeof(int64_t));
        if ((*F_handle)->parent == NULL)
        {
            SPEX_factorization_free(F_handle);
            return SPEX_OUT_OF_MEMORY;
        }
        memcpy((*F_handle)->parent, S->parent, A->n * sizeof(int64_t));
        */
    }
    else // TODO QR
    {
        // set the kind for factorization
        (*F_handle)->kind = SPEX_QR_FACTORIZATION;
    }
    
    return SPEX_OK;
}
