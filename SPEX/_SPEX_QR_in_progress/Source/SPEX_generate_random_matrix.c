//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_generate_random_matrix.c: Generate a random dense matrix
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


# include "spex_qr_internal.h"

SPEX_info SPEX_generate_random_matrix
(
    SPEX_matrix *A_handle, // Matrix to be created. Null on input
    int64_t m,              // Rows of the matrix
    int64_t n,              // Columns of the matrix
    unsigned int seed,      // Random number seed
    int64_t lower,          // Lower bound for numbers to be generated
    int64_t upper           // Upper bound for numbers to be generated
)
{
    // Input checks
    ASSERT(n >= 0);
    ASSERT(m >= 0);
    ASSERT(lower < upper);
    
    SPEX_info info;
    SPEX_matrix A = NULL; 
    // A is a m*n triplet matrix whose entries are FP64 Note that the first
    // boolean parameter says that the matrix is not shallow, so that A->i,
    // A->j, and A->x are calloc'd. The second boolean parameter is meaningless
    // for FP64 matrices, but it tells SPEX to allocate the values of A->x
    // for the mpz_t, mpq_t, and mpfr_t entries
    SPEX_CHECK(SPEX_matrix_allocate(&A, SPEX_TRIPLET, SPEX_FP64, m, n, m*n,
        false, true, NULL));
    
    int64_t counter = 0;
    // Start the random number generator with seed
    srand(seed);
    for (int64_t k = 0; k < m; k++)
    {
        for (int64_t p = 0; p < n; p++)
        {
            A->i[counter] = k;
            A->j[counter] = p;
            A->x.fp64[counter] = rand() % (upper-lower) + lower;
            counter+=1;
        }
    }
    (*A_handle) = A;
    return SPEX_OK;
}
