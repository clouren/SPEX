//------------------------------------------------------------------------------
// SPEX_Util/SPEX_matrix_nnz: find # of entries in a matrix
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#pragma GCC diagnostic ignored "-Wunused-variable"
#include "spex_util_internal.h"

int64_t SPEX_matrix_nnz     // return # of entries in A, or -1 on error
(
    const SPEX_matrix *A,      // matrix to query
    const SPEX_options *option // command options, currently unused
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!spex_initialized ( )) return (-1) ;

    if (A == NULL)
    {
        return (-1) ;
    }

    //--------------------------------------------------------------------------
    // find nnz (A)
    //--------------------------------------------------------------------------

    // In all three cases, SPEX_matrix_nnz(A,option) is <= A->nzmax.

    switch (A->kind)
    {
        case SPEX_CSC:
            // CSC matrices:  nnz(A) is given by Ap[n].  A->nz is ignored.
            return ((A->p == NULL || A->n < 0) ? (-1) : A->p [A->n]) ;
        case SPEX_TRIPLET:
            // triplet matrices:  nnz(A) is given by A->nz.
            return (A->nz) ;
        case SPEX_DENSE:
            // dense matrices: nnz(A) is always m*n.  A->nz is ignored.
            return ((A->m < 0 || A->n < 0)? (-1) : (A->m * A->n)) ;
        default:
            return (-1) ;
    }
}

