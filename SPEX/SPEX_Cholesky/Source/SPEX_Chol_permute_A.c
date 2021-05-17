//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_permute_A: Symmetric permutation of matrix A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

/* Purpose: Permute the matrix A and return A2 = PAP */
SPEX_info SPEX_Chol_permute_A
(
    SPEX_matrix** A2_handle, // Output permuted matrix
    SPEX_matrix* A,          // Initial input matrix
    SPEX_Chol_analysis* S    //Symbolic analysis struct that contains column 
                            //and inverse row permutations
)
{
    SPEX_info info;
    // Check inputs
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!A2_handle || !S)
        return SPEX_INCORRECT_INPUT;

    //Create pinv, the inverse row permutation
    int k, index;
    int n = A->n;
    int64_t* pinv = NULL;

    pinv = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv[index] = k;
    }
    S->pinv=pinv;
   

    // Allocate memory for A2 which is a permuted copy of A
    int64_t nz = 0, j;
    SPEX_matrix* A2 = NULL;
    SPEX_matrix_allocate(&A2, SPEX_CSC, SPEX_MPZ, n, n, A->p[A->n], false, true, NULL);


    // Set A2 scale
    SPEX_CHECK(SPEX_mpq_set(A2->scale, A->scale));
    // Populate the entries in A2
    for (int64_t k = 0 ; k < n ; k++)
    {
        A2->p [k] = nz ;                       // column k of A2 is column q[k] of A 
        j = S->q [k];
        for (int64_t t = A->p [j] ; t < A->p [j+1] ; t++)
        {
            SPEX_CHECK(SPEX_mpz_set(A2->x.mpz[nz], A->x.mpz[t]));  // row i of A is row pinv[i] of A2 
            A2->i [nz++] = pinv [A->i [t]];
        }
    }
    A2->p [n] = nz ;                       // finalize the last column of A2 
    (*A2_handle) = A2;                     // Set A2_handle
    return SPEX_OK;
}
