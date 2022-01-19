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


#include "spex_update_internal.h"

#define SPEX_FREE_ALL   \
    SPEX_factorization_free(&F, option);

SPEX_info SPEX_Update_factorization_convert
(
    SPEX_factorization **F_out,// The output factorization with same
                            // factorization kind as F_in
    const SPEX_factorization *F_in, // The factorization to be converted
    const bool updabtable, // true if wish to obtain updatable F.
    const SPEX_options* option // Command options
)  
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    if (!spex_initialized()) {return SPEX_PANIC;}
    
    if (!F_out || !F_in ||
        (F_in->kind != SPEX_LU_FACTORIZATION &&
         F_in->kind != SPEX_CHOLESKY_FACTORIZATION))
    {
        return SPEX_INCORRECT_INPUT;
    }

    (*F_out) = NULL;

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------
    SPEX_info info;
    SPEX_factorization *F = NULL;
    int64_t i, n = F_in->L->n;

    //--------------------------------------------------------------------------
    // allocate memory space for the factorization
    //--------------------------------------------------------------------------
    F = (SPEX_factorization*) SPEX_calloc(1, sizeof(SPEX_factorization));
    if (F == NULL)
    {
        return SPEX_OUT_OF_MEMORY;
    }
    // set factorization kind
    F->kind = F_in->kind;
    // Allocate and set scale_for_A
    SPEX_CHECK(SPEX_mpq_init(F->scale_for_A));
    SPEX_CHECK(SPEX_mpq_set (F->scale_for_A, F_in->scale_for_A));

    //--------------------------------------------------------------------------
    // row permutations
    //--------------------------------------------------------------------------
    F->P_perm    = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
    F->Pinv_perm = (int64_t*) SPEX_malloc (n * sizeof(int64_t));

    if (!(F->P_perm) || !(F->Pinv_perm))
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    memcpy(F->P_perm,    F_in->P_perm,     n * sizeof(int64_t));
    memcpy(F->Pinv_perm, F_in->Pinv_perm,  n * sizeof(int64_t));

    //--------------------------------------------------------------------------
    // obtain column permutations for LU factorization
    //--------------------------------------------------------------------------
    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        F->Q_perm    = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
        if (!(F->Q_perm))
        {
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }
        memcpy(F->Q_perm,    F_in->Q_perm,     n * sizeof(int64_t));

        // obtain Qinv_perm only for updatable LU factorization
        if (updatable)
        {
            // generate Qinv_perm for F
            F->Qinv_perm = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
            if (!(F->Qinv_perm))
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }
            for (i = 0; i < n; i++)
            {
                F->Qinv_perm[F->Q_perm[i]] = i;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Create rhos, a global dense mpz_t matrix of dimension n*1
    //--------------------------------------------------------------------------
    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, false, option));
    // copy rhos
    for (i = 0; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpz_set(SPEX_1D(F->rhos,    i, mpz),
                                SPEX_1D(F_in->rhos, i, mpz));
    }

    //--------------------------------------------------------------------------
    // obtain the desired format of matrices
    //--------------------------------------------------------------------------
    if (updatable)
    {
        SPEX_CHECK(SPEX_matrix_copy(&(F->L), SPEX_DYNAMIC_CSC, SPEX_MPZ,
            F_in->L, option));
        SPEX_CHECK(SPEX_Update_permute_row(F->L, F->P_perm, option));
        SPEX_CHECK(SPEX_Update_matrix_canonicalize(F->L, F->P_perm, option));
        if (F->kind == SPEX_LU_FACTORIZATION)
        {
            SPEX_matrix *UT = NULL;
            SPEX_CHECK(SPEX_transpose(&UT, F->U, option));
            SPEX_CHECK(SPEX_matrix_copy(&(F->U), SPEX_DYNAMIC_CSC, SPEX_MPZ,
                UT, option));
            SPEX_CHECK(SPEX_matrix_free(&UT, option));
            SPEX_CHECK(SPEX_Update_permute_row(F->U, F->Q_perm, option));
            SPEX_CHECK(SPEX_Update_matrix_canonicalize(F->U, F->Q_perm,
                option));
        }
    }
    else
    {
        SPEX_CHECK(SPEX_matrix_copy(&(F->L), SPEX_CSC, SPEX_MPZ,
            F_in->L, option));
        SPEX_CHECK(SPEX_Update_permute_row(F->L, F->Pinv_perm, option));
        if (F->kind == SPEX_LU_FACTORIZATION)
        {
            SPEX_matrix *U_CSC = NULL;
            SPEX_CHECK(SPEX_matrix_copy(&U_CSC, SPEX_CSC, SPEX_MPZ,
                F->U, option));
            SPEX_CHECK(SPEX_transpose(&(F->U), U_CSC, option));
            SPEX_CHECK(SPEX_matrix_free(&U_CSC, option));
            SPEX_CHECK(SPEX_Update_permute_row(F->U, F->Qinv_perm, option));
        }
    }

    (*F_out) = F;
    return SPEX_OK;
}
