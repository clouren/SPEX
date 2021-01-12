//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_symbolic_add: Symbolically add two sparse matrices
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE      \
    SPEX_FREE(x);           \
    SPEX_FREE(w);           \

#include "SPEX_WAMF.h"


/* Purpose: C = A + B. This function gives a rough upper bound in 
 * the number of bits in A+B. 
 * 
 * On success, the matrix C is returned.
 * 
 */

SPEX_matrix* SPEX_WAMF_symbolic_add
(
    SPEX_matrix *A,      // Left matrix
    SPEX_matrix *B      // Right matrix
)
{
    if (!A || !A->p || !A->i || !A->x.fp64 || !B || !B->p || !B->i || !B->x.fp64)
        return NULL;
 
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_FP64);
    
    ASSERT(B->kind == SPEX_CSC);
    ASSERT(B->type == SPEX_FP64);
    
    /* Allocate memory */
    int64_t p, j, nz = 0, m, n;
    m = A->m ; n = B->n ;
    SPEX_matrix *C = NULL;
    SPEX_matrix_allocate(&C, SPEX_CSC, SPEX_FP64, m, n, A->nzmax+B->nzmax,
        false, true, NULL);
        
    double* x = SPEX_calloc(m, sizeof(double));
    int64_t* w = SPEX_calloc(m, sizeof(int64_t));
    if (!C || !x || !w)
    {
        FREE_WORKSPACE;
        SPEX_matrix_free(&C, NULL);
        return NULL;
    }
        
    /* Add C = A+B */
    for (j = 0 ; j < n ; j++)
    {
        C->p [j] = nz ;                                         /* column j of C starts here */
        nz = SPEX_WAMF_scatter (A, j, 1, w, x, j+1, C, nz) ;       /* A(:,j)*/
        nz = SPEX_WAMF_scatter (B, j, 1, w, x, j+1, C, nz) ;        /* B(:,j) */
        for (p = C->p [j] ; p < nz ; p++) C->x.fp64[p] = x[C->i[p]];
    }
    SPEX_options* option;
    C->p [n] = nz ; C->nz = nz;                                 /* finalize the last column of C */
    /* Free memory */
    FREE_WORKSPACE;
    // Collapse the matrix C
    bool ok;
    C->i = (int64_t*) SPEX_realloc( A->p[n], A->nzmax, sizeof(int64_t), C->i, &ok);
    C->x.fp64 = (double*) SPEX_realloc( A->p[n], A->nzmax, sizeof(double), C->x.fp64, &ok);
    C->nzmax = C->p[C->n];
    return C;
}
