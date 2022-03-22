//------------------------------------------------------------------------------
// SPEX_Update/spex_update_factorization_convert.c: convert between updatable
// and non-updatable factorization.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// spex_update_factorization_convert is called either to obtain the updatable
// factorization from non-updatable factorization with matrices in any kind or
// type, or to obtain non-updatable factorization with SPEX_CSC matrices with
// MPZ entries.
// To obtain the updatable format, this function will obtain the
// transpose of U (for LU factorization), permute rows of L and UT, and make
// sure each column of the matrices have cooresponding pivot as the first
// entry. To otain the non-updatable format, this function will transpose UT
// (for LU factorization) and permute rows and L and U.

#define SPEX_FREE_ALL             \
    SPEX_FREE(Qinv);              \
    SPEX_matrix_free(&L, option); \
    SPEX_matrix_free(&U, option); \
    SPEX_matrix_free(&Mat, option);


#include "spex_update_internal.h"

SPEX_info spex_update_factorization_convert
(
    SPEX_factorization *F, // The factorization to be converted
    const bool updatable, // true if wish to obtain updatable F.
    const SPEX_options* option // Command options
)  
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    if (!spex_initialized()) {return SPEX_PANIC;}
    
    if (!F ||
        (F->kind != SPEX_LU_FACTORIZATION &&
         F->kind != SPEX_CHOLESKY_FACTORIZATION))
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------
    SPEX_info info;
    int64_t i, n = F->L->n;
    int64_t *Qinv = NULL;
    SPEX_matrix *L = NULL, *U = NULL, *Mat = NULL;

    //--------------------------------------------------------------------------
    // obtain Qinv_perm for updatable LU factorization
    //--------------------------------------------------------------------------
    if (F->kind == SPEX_LU_FACTORIZATION && updatable)
    {
        Qinv = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
        if (!Qinv)
        {
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }
        for (i = 0; i < n; i++)
        {
            Qinv[F->Q_perm[i]] = i;
        }
        SPEX_FREE(F->Qinv_perm); // clear whatever is there just in case
        F->Qinv_perm = Qinv;
        Qinv = NULL;
    }

    //--------------------------------------------------------------------------
    // obtain the desired format of L and/or U
    //--------------------------------------------------------------------------
    if (updatable)
    {
        SPEX_CHECK(SPEX_matrix_copy(&L, SPEX_DYNAMIC_CSC, SPEX_MPZ,
            F->L, option));
        SPEX_CHECK(spex_update_permute_row(L, F->P_perm, option));
        SPEX_CHECK(spex_update_matrix_canonicalize(L, F->P_perm, option));
        // replace the original L
        SPEX_CHECK(SPEX_matrix_free(&(F->L), option));
        F->L = L;
        L = NULL;

        if (F->kind == SPEX_LU_FACTORIZATION)
        {
            // TODO do the shallow copy of all mpz values
            SPEX_CHECK(SPEX_transpose(&Mat, F->U, option));
            SPEX_CHECK(SPEX_matrix_copy(&U, SPEX_DYNAMIC_CSC, SPEX_MPZ,
                Mat, option));
            SPEX_CHECK(spex_update_permute_row(U, F->Q_perm, option));
            SPEX_CHECK(spex_update_matrix_canonicalize(U, F->Q_perm,
                option));
            // replace the original U
            SPEX_CHECK(SPEX_matrix_free(&(F->U), option));
            F->U = U;
            U = NULL;
        }
    }
    else
    {
        SPEX_CHECK(SPEX_matrix_copy(&L, SPEX_CSC, SPEX_MPZ, F->L, option));
        SPEX_CHECK(spex_update_permute_row(L, F->Pinv_perm, option));
        // replace the original L
        SPEX_CHECK(SPEX_matrix_free(&(F->L), option));
        F->L = L;
        L = NULL;

        if (F->kind == SPEX_LU_FACTORIZATION)
        {
            SPEX_CHECK(SPEX_matrix_copy(&Mat, SPEX_CSC, SPEX_MPZ,
                F->U, option));
            SPEX_CHECK(SPEX_transpose(&U, Mat, option));
            SPEX_CHECK(spex_update_permute_row(U, F->Qinv_perm, option));
            // replace the original U
            SPEX_CHECK(SPEX_matrix_free(&(F->U), option));
            F->U = U;
            U = NULL;
        }
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
