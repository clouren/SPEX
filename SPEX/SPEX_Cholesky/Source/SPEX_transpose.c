//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_transpose: Transpose a matrix
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE  \
    SPEX_FREE(w);       \

#include "spex_chol_internal.h"
    
//TODO: Move to SPEX_Util
/* Purpose: This function sets C = A' 
 * C_handle is NULL on input. On output, C_handle contains a pointer to A'
 */
SPEX_info SPEX_transpose
(
    SPEX_matrix **C_handle,     // C = A'
    SPEX_matrix *A              // Matrix to be transposed
)
{
    SPEX_info ok;
    // Check input
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!C_handle) 
        return SPEX_INCORRECT_INPUT;
    
    // Declare workspace and C
    int64_t* w = NULL;
    SPEX_matrix* C = NULL;
    OK(SPEX_matrix_allocate(&C, SPEX_CSC, SPEX_MPZ, A->n, A->m, A->p[A->n], false, true, NULL));
    int64_t p, q, j, n, m;
    m = A->m ; n = A->n ; 
    w = (int64_t*) SPEX_malloc(m* sizeof(int64_t));
    if (!w)
    {
        SPEX_matrix_free(&C, NULL);
        return SPEX_OUT_OF_MEMORY;
    }
    for (p = 0; p < m; p++) w[p] = 0;
    for (p = 0 ; p < A->p [n] ; p++) w [A->i [p]]++ ;       // row counts 
    SPEX_cumsum (C->p, w, m) ;                           // row pointers
    for (j = 0 ; j < n ; j++)
    {
        for (p = A->p [j] ; p < A->p [j+1] ; p++)
        {
            q = w [A->i [p]]++;
            C->i [q] = j ;                 // place A(i,j) as entry C(j,i) 
            OK(SPEX_mpz_set(C->x.mpz[q], A->x.mpz[p]));
        }
    }
    C->p[m] = A->p[n];
    (*C_handle) = C;
    FREE_WORKSPACE;
    return SPEX_OK;
}
