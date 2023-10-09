//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_transpose: Transpose a CSC matrix
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORK       \
    SPEX_FREE(w);

#define SPEX_FREE_ALL        \
    SPEX_FREE_WORK;          \
    SPEX_matrix_free(&C, option);

#include "spex_util_internal.h"

/* Purpose: This function sets C = A', where A must be a SPEX_CSC matrix
 * C_handle is NULL on input. On output, C_handle contains a pointer to A'
 */

SPEX_info spex_qr_transpose
(
    SPEX_matrix *C_handle,      // C = A'
    SPEX_matrix A,              // Matrix to be transposed
    const SPEX_options option
)
{

    SPEX_info info;
    if (!spex_initialized ( )) return (SPEX_PANIC);
    // Check input
    SPEX_REQUIRE_KIND (A, SPEX_CSC);
    if (!C_handle)       { return SPEX_INCORRECT_INPUT;}

    // Declare workspace and C
    int64_t *w = NULL;
    SPEX_matrix C = NULL;
    int64_t nz;                            // Number of nonzeros in A
    int64_t p, q, j, n, m;
    info = SPEX_matrix_nnz(&nz, A, option);
    if (info != SPEX_OK) {return info;}
    m = A->m ; n = A->n ;
    ASSERT( m >= 0);
    ASSERT( n >= 0);

    // C is also CSC and its type is the same as A
    SPEX_CHECK(SPEX_matrix_allocate(&C, SPEX_CSC, A->type, n, m, nz,
        false, true, option));

    // Declare workspace
    w = (int64_t*) SPEX_calloc(m, sizeof(int64_t));
    if (!w)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    // Compute row counts
    for (p = 0 ; p < nz ; p++)
    {
        w [A->i [p]]++ ;
    }

    // Compute row pointers
    spex_cumsum (C->p, w, m);

    for (j = 0 ; j < n ; j++)
    {
        for (p = A->p [j] ; p < A->p [j+1] ; p++)
        {
            q = w [A->i [p]]++;
            C->i [q] = j ;                 
        }
    }
    //SPEX_MPQ_SET(C->scale, A->scale);

    (*C_handle) = C;
    SPEX_FREE_WORK;
    return SPEX_OK;
}
