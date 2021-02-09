//------------------------------------------------------------------------------
// SPEX_Util/SPEX_transpose: Transpose a matrix
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// TODO Delete me for release

#define FREE_WORKSPACE  \
    SPEX_FREE(w);       \

#include "spex_util_internal.h"
    
/* Purpose: This function sets C = A' 
 * C_handle is NULL on input. On output, C_handle contains a pointer to A'
 */
SPEX_info SPEX_transpose
(
    SPEX_matrix **C_handle,     // C = A'
    SPEX_matrix *A              // Matrix to be transposed
)
{
    SPEX_info info;
    // Check input
    ASSERT(A->kind == SPEX_CSC);
    if (!C_handle) 
        return SPEX_INCORRECT_INPUT;
    
    // Declare workspace and C
    int64_t* w = NULL;
    SPEX_matrix* C = NULL;
    // C is also CSC and its type is the same as A
    SPEX_CHECK(SPEX_matrix_allocate(&C, SPEX_CSC, A->type, A->n, A->m, A->p[A->n], false, true, NULL));
    int64_t p, q, j, n, m;
    m = A->m ; n = A->n ; 
    
    ASSERT( m >= 0);
    ASSERT( n >= 0);
    // Declare workspace
    w = (int64_t*) SPEX_calloc(m, sizeof(int64_t));
    if (!w)
    {
        SPEX_matrix_free(&C, NULL);
        return SPEX_OUT_OF_MEMORY;
    }
    // Compute row counts
    for (p = 0 ; p < A->p [n] ; p++) w [A->i [p]]++ ;
    // Compute row pointers
    SPEX_cumsum (C->p, w, m) ;
    // Populate C
    for (j = 0 ; j < n ; j++)
    {
        for (p = A->p [j] ; p < A->p [j+1] ; p++)
        {
            q = w [A->i [p]]++;
            C->i [q] = j ;                 // place A(i,j) as entry C(j,i) 
            // Place the values. The way values are set are based on the type of A
            // TODO is there a better way to do this?
            if (A->type == SPEX_MPZ)
            {
                SPEX_CHECK(SPEX_mpz_set(C->x.mpz[q], A->x.mpz[p]));
            }
            else if (A->type == SPEX_MPQ)
            {
                SPEX_CHECK(SPEX_mpq_set(C->x.mpq[q], A->x.mpq[p]));
            }
            else if (A->type == SPEX_MPFR)
            {
                SPEX_CHECK(SPEX_mpfr_set(C->x.mpfr[q], A->x.mpfr[p], SPEX_DEFAULT_MPFR_ROUND));
            }
            else if (A->type == SPEX_INT64)
            {
                C->x.int64[q] = A->x.int64[p];
            }
            else
                C->x.fp64[q] = A->x.fp64[p];
        }
    }
    C->p[m] = A->p[n];
    (*C_handle) = C;
    FREE_WORKSPACE;
    return SPEX_OK;
}
