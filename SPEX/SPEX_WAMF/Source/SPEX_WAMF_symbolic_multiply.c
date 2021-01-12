//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_symbolic_multiply: Symbolically multiply two WAMF_sparse matrices
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: C = A*B. This function gives a rough upper bound in 
 * the number of bits in A*B. Note that A and B are not explicitly multiplied since
 * we are dealing with symbolic bit operations
 * 
 * On success, the matrix C is returned.
 * 
 */

SPEX_matrix* SPEX_WAMF_symbolic_multiply
(
    SPEX_matrix *A,  // Left matrix
    SPEX_matrix *B   // Right matrix
)
{
    if (!A || !A->x.fp64 || !A->i || !A->p || !B || !B->p || !B->x.fp64 || !B->i)
        return NULL;
    
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_FP64);
    
    ASSERT(B->kind == SPEX_CSC);
    ASSERT(B->type == SPEX_FP64);
    
    /* Allocate memory */
    int64_t p, j, nz = 0, m, n ;
    m = A->m; n = B->n;
    double *x = SPEX_calloc(m, sizeof(double));
    int64_t* w = SPEX_calloc(m, sizeof(int64_t));
    SPEX_matrix* C = NULL;
    SPEX_matrix_allocate(&C, SPEX_CSC, SPEX_FP64, m, n, A->nzmax+B->nzmax,
        false, true, NULL);
    
    
    if (!x || !w || !C)
    {
        SPEX_FREE(x); SPEX_FREE(w); SPEX_matrix_free(&C, NULL);
        return NULL;
    }
    
    /* Multiply A*B */
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax )         // If need more nonzeros
        {
            SPEX_WAMF_sparse_realloc(C);
            C->nz = nz;
        }
        C->p [j] = nz ;                   /* column j of C starts here */
        for (p = B->p [j] ; p < B->p [j+1] ; p++)
        {
            nz = SPEX_WAMF_scatter (A, B->i [p], B->x.fp64 [p], w, x, j+1, C, nz) ;
        }
        for (p = C->p [j] ; p < nz ; p++) C->x.fp64[p] = x[C->i[p]];
        C->nz = nz;
    }
    C->p [n] = nz ; C->nz = nz;
    SPEX_FREE(w); SPEX_FREE(x);
    // Remove extra space from C
    bool ok;
    C->i = (int64_t*) SPEX_realloc( A->p[n], A->nzmax, sizeof(int64_t), C->i, &ok);
    C->x.fp64 = (double*) SPEX_realloc( A->p[n], A->nzmax, sizeof(double), C->x.fp64, &ok);
    C->nzmax = C->p[C->n];
    return C ;  
}
