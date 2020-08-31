//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_permute_A: Symmetric permutation of matrix A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_Chol.h"

/* Purpose: Permute the matrix A and return A2 = PAP */
SPEX_info SPEX_Chol_permute_A
(
    SPEX_matrix **A2_handle, // Output permuted matrix
    SPEX_matrix* A,          // Initial input matrix
    int64_t* pinv,           // Inverse row permutation
    SPEX_LU_analysis* S      // Column permutation
)
{
    SPEX_info ok;
    // Check inputs
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!A2_handle || !pinv || !S)
        return SPEX_INCORRECT_INPUT;
    
    // Allocate memory for A2 which is a permuted copy of A
    int64_t nz = 0, j, n = A->n;
    SPEX_matrix* A2 = NULL;
    SPEX_matrix_allocate(&A2, SPEX_CSC, SPEX_MPZ, n, n, A->p[A->n], false, true, NULL);

    // Set A2 scale
    OK(SPEX_mpq_set(A2->scale, A->scale));
    // Populate the entries in A2
    for (int64_t k = 0 ; k < n ; k++)
    {
        A2->p [k] = nz ;                       // column k of A2 is column q[k] of A 
        j = S->q [k];
        for (int64_t t = A->p [j] ; t < A->p [j+1] ; t++)
        {
            OK(SPEX_mpz_set(A2->x.mpz[nz], A->x.mpz[t]));  // row i of A is row pinv[i] of A2 
            A2->i [nz++] = pinv [A->i [t]];
        }
    }
    A2->p [n] = nz ;                       // finalize the last column of A2 
    (*A2_handle) = A2;                     // Set A2_handle
    return SPEX_OK;
}
