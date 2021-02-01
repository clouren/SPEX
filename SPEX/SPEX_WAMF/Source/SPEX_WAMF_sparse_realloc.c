//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_sparse_realloc: Double size of a WAMF_sparse matrix
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"


/* Purpose: This function doubles the size of a WAMF_sparse matrix. 
 */

SPEX_info SPEX_WAMF_sparse_realloc
(
    SPEX_matrix* A
)
{
    if (!A || !A->p || !A->i){return SPEX_INCORRECT_INPUT;}
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_FP64);
    int64_t nzmax = A->nzmax;
    // Double size of A->x and A->i
    bool ok;
    A->i = (int64_t*) SPEX_realloc(2*nzmax, nzmax, sizeof(int64_t), A->i, &ok);
    A->x.fp64 = (double*) SPEX_realloc(2*nzmax, nzmax, sizeof(double), A->x.fp64, &ok);
    A->nzmax = nzmax*2;
    
    if (!A->i || !A->x.fp64) {return SPEX_OUT_OF_MEMORY;}

    return SPEX_OK;
}
