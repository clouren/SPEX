//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_scatter: Sparse scatter x = x + A(:,j)
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse 
 * This function gives a rough upper bound of the number of bits present in x
 */

int64_t SPEX_WAMF_scatter
(
    SPEX_matrix *A, // Input matrix
    int64_t j,          // Column of A
    double beta,    //
    int64_t *w,         // workspace vector
    double *x,      // x = x + beta A(:,j)
    int64_t mark,       // location of C
    SPEX_matrix *C, // Set in C
    int64_t nz          // Number of nonzeros
)
{
    int64_t i, p;
    
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_FP64);
    
    ASSERT(C->kind == SPEX_CSC);
    ASSERT(C->type == SPEX_FP64);
    
    for (p = A->p [j] ; p < A->p [j+1] ; p++)
    {
        i = A->i [p] ;                          /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            C->i [nz++] = i ;                   /* add i to pattern of C(:,j) */
            x[i] = beta + A->x.fp64[p];        // Number of bits in x[i] <= bits in beta + bits in A->x[p]
        }
        else
            x[i] = SPEX_MAX( x[i], A->x.fp64[p]) + 1;     // Number of bits in x[i] <= max( bits in x[i], bits in beta*Ax->[p]) + 1
    }
    return nz ;
}
