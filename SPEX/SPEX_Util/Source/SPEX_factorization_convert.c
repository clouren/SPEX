//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_factorization_convert.c: convert between updatable
// and non-updatable factorization.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// spex_update_factorization_convert is called to make conversion between the
// updatable factorization (MPZ & DYNAMIC_CSC) and non-updatable MPZ & CSC
// factorization. That is, when the input factorization is updatable with MPZ
// entries and DYNAMIC_CSC kind, the returned factorization will be
// non-updatable with MPZ entries and CSC kind. While if the input
// factorization is in non-updatable format with MPZ entries and CSC kind, the
// returned factorization will be updatable with MPZ entries and DYNAMIC_CSC
// kind. If the matrices of the factorization are in kinds other than CSC or
// DYNAMIC_CSC, or in types other than MPZ, the function return INCORRECT_INPUT
// error.
// To obtain the updatable format, this function will obtain the
// transpose of U (for LU factorization), permute rows of L and UT, and make
// sure each column of the matrices have cooresponding pivot as the first
// entry. To otain the non-updatable format, this function will transpose UT
// (for LU factorization) and permute rows and L and U.
// If the function return unsuccessfully, the returned factorization should be
// considered as undefined.
//
// NOTE: if F->updatable == false upon input, F->L (and F->U if exists) must be
// CSC MPZ matrix, otherwise, SPEX_INCORRECT_INPUT will be returned. Likewise,
// if F->updatable == true upon input, F->L (and F->U if exists) must be
// dynamic_CSC MPZ matrix. In addition, both F->L and F->U (if exists) must not
// be shallow matrices. All SPEX functions output either of these two formats
// and non-shallow. Therefore, these input requirements can be met easily if
// users do not try to modify any individual component of F.  The conversion is
// done in place and F->updatable will be set to its complement upon output. In
// case of any error, the returned factorization should be considered as
// undefined.



#include "spex_util_internal.h"

SPEX_info SPEX_factorization_convert
(
    SPEX_factorization *F, // The factorization to be converted
    const SPEX_options* option // Command options
)  
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    if (!spex_initialized()) {return SPEX_PANIC;}
    if (!F || !(F->P_perm) || !(F->Pinv_perm) || F->L->type != SPEX_MPZ ||
        (F->L->kind != SPEX_CSC && F->L->kind != SPEX_DYNAMIC_CSC) ||
        (F->kind != SPEX_LU_FACTORIZATION &&
         F->kind != SPEX_CHOLESKY_FACTORIZATION) ||
        (F->kind == SPEX_LU_FACTORIZATION && (!(F->Q_perm) ||
         (F->updatable && !(F->Qinv_perm)) || F->U->type != SPEX_MPZ ||
         (F->U->kind != SPEX_CSC && F->U->kind != SPEX_DYNAMIC_CSC))))
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------
    SPEX_info info;
    int64_t i, n = F->L->n;

    // update the updatable flag
    F->updatable = !(F->updatable);

    //--------------------------------------------------------------------------
    // obtain/update Qinv_perm for updatable LU factorization
    //--------------------------------------------------------------------------
    if (F->kind == SPEX_LU_FACTORIZATION && F->updatable)
    {
        // FIXME is this ok?
        // Although Qinv_perm is NULL when F is created by factorizing matrix,
        // where F is initially not updatable, Qinv_perm is then created when
        // converted to updatable format and not deleted even when converted
        // back to non-updatable from updatable format.

        // Allocate space for Qinv_perm if it was NULL. Otherwise, Qinv_perm
        // is re-usable since the size is never changed.
        if (!(F->Qinv_perm))
        {
            F->Qinv_perm = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
            if (!(F->Qinv_perm))
            {
                return SPEX_OUT_OF_MEMORY;
            }
        }
        for (i = 0; i < n; i++)
        {
            F->Qinv_perm[F->Q_perm[i]] = i;
        }
    }

    //--------------------------------------------------------------------------
    // obtain the desired format of L and/or U
    //--------------------------------------------------------------------------
    // convert matrix L
    SPEX_CHECK(spex_matrix_convert(F, true, option));

    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        // convert matrix U
        SPEX_CHECK(spex_matrix_convert(F, false, option));
    }

    return SPEX_OK;
}
