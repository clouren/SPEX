//------------------------------------------------------------------------------
// SPEX_Util/SPEX_symbolic_analyze: symbolic ordering and analysis for specified
// kind
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* SPEX_symbolic_analyze performs the symbolic ordering and analysis for
 * specified factorization kind. The pre-ordering method is specified in the
 * option.
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

// SPEX_symbolic_analyze creates the SPEX_symbolic_analysis object S.  Use
// SPEX_symbolic_analysis_free to delete it.

#include "spex_util_internal.h"

SPEX_info SPEX_symbolic_analyze
(
    SPEX_symbolic_analysis **S_handle, // symbolic analysis
    const SPEX_matrix *A, // Input matrix
    const SPEX_symbolic_analysis_kind kind, // LU, Cholesky or QR analysis
    const SPEX_options *option  // Control parameters
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    if (!spex_initialized()) return SPEX_PANIC;

    // A can have any data type, but must be in sparse CSC format
    SPEX_REQUIRE_KIND (A, SPEX_CSC) ;
    
    if (!S_handle || A->n != A->m ||
        (kind != SPEX_LU_SYMBOLIC_ANALYSIS &&
         kind != SPEX_CHOLESKY_SYMBOLIC_ANALYSIS &&
         kind != SPEX_QR_SYMBOLIC_ANALYSIS))
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*S_handle) = NULL;
    SPEX_info info;

    //--------------------------------------------------------------------------
    // perform symbolic analysis for specified kind
    //--------------------------------------------------------------------------
    if (kind == SPEX_LU_SYMBOLIC_ANALYSIS)
    {
        // perform symbolic analysis for LU factorization
        info = SPEX_LU_analyze(S_handle, A, option);
        if (info != SPEX_OK)       {  return info;  }
    }
    else if (kind == SPEX_CHOLESKY_SYMBOLIC_ANALYSIS)
    {
        // perform symbolic analysis for Cholesky factorization
        // TODO this is in SPEX_Cholesky
        // TODO remove input checking that has been done here
        // TODO set SPEX_factorization kind
        info = SPEX_Chol_preorder(S_handle, A, option);
        if (info != SPEX_OK)      {  return info;   }
    }
    else // if (kind == SPEX_QR_ANALYSIS)
    {
        // TODO
    }

    return SPEX_OK;
}
