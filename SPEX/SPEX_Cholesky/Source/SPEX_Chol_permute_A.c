//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_permute_A: Symmetric permutation of matrix A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

/* Purpose: Permute the matrix A and return PAP = PAP 
 * Input arguments:
 * 
 * PAP_handle:   The user's permuted input matrix. 
 * 
 * A:           The user's input matrix
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. 
 *              Contains column and inverse row permutations
 * 
 */
SPEX_info SPEX_Chol_permute_A
(
    //Output
    SPEX_matrix** PAP_handle,    // Output permuted matrix
    //Input
    SPEX_matrix* A,             // Input matrix
    SPEX_Chol_analysis* S      //Symbolic analysis struct that contains column 
                             //and inverse row permutations
)
{
    SPEX_info info;
    // Check inputs
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!PAP_handle || !S)
        return SPEX_INCORRECT_INPUT;

    //Create pinv, the inverse row permutation
    int k, index;
    int n = A->n;
    int64_t* pinv = NULL;

    pinv = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    for (k = 0; k < n; k++)
    {
        index = S->p[k];
        pinv[index] = k;
    }
    S->pinv=pinv;
    

    // Allocate memory for PAP which is a permuted copy of A
    int64_t nz = 0, j;
    SPEX_matrix* PAP = NULL;
    SPEX_matrix_allocate(&PAP, SPEX_CSC, SPEX_MPZ, n, n, A->p[A->n], false, true, NULL);


    // Set PAP scale
    SPEX_CHECK(SPEX_mpq_set(PAP->scale, A->scale));
    // Populate the entries in PAP
    for (int64_t k = 0 ; k < n ; k++)
    {
        PAP->p [k] = nz ;                       // column k of PAP is column p[k] of A 
        j = S->p [k];
        for (int64_t t = A->p [j] ; t < A->p [j+1] ; t++)
        {
            SPEX_CHECK(SPEX_mpz_set(PAP->x.mpz[nz], A->x.mpz[t]));  // row i of A is row pinv[i] of PAP 
            PAP->i [nz++] = pinv [A->i [t]];
        }
    }
    PAP->p [n] = nz ;                       // finalize the last column of PAP
    (*PAP_handle) = PAP;                     // Set PAP_handle
    return SPEX_OK;
}
