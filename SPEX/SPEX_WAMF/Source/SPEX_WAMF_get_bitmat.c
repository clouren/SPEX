//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_get_bitmat: Get a matrix of bit-lengths
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: Given a matrix A which is stored as a sparse CSC matrix with GMP entries, 
 * convert it to a matrix of bit-lengths. That is the output matrix, B, has the same nonzero
 * pattern as A, but the entries are B(i,j) = bit-length( A(i,j)).
 * 
 * On success, this function returns the new matrix B. On failure, this function returns
 * a NULL pointer. 
 */

SPEX_matrix* SPEX_WAMF_get_bitmat
(
    SPEX_matrix *A
)
{
    if ( !A || !A->x.mpz || !A->i || !A->p ) {return NULL;}
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);
    
    // Create matrix B
    SPEX_matrix* B = NULL;
    SPEX_matrix_allocate(&B, SPEX_CSC, SPEX_FP64, A->m, A->n, A->nzmax,
        false, true, NULL);
    
    // Create the values of B    
    if (!B) return NULL;
    for (int64_t i = 0; i <= A->n; i++)             // Set B->p
        B->p[i] = A->p[i];
    for (int64_t i = 0; i < A->p[A->n]; i++)             // Set B->i and B->x
    {
        B->i[i] = A->i[i];                      // Set B->i
        size_t size;
        SPEX_mpz_sizeinbase( &size, A->x.mpz[i], 2);    // Number of bits in A->x[i]
        B->x.fp64[i] = (double) size;
    }
    return B;
}
