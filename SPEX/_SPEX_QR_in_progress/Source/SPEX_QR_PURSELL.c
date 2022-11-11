//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_QR_PURSELL.c: Pursell Algorithm
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs REF QR via the Pursell Algorithm on A'*A.
 */

# include "spex_qr_internal.h"

//* Perform the IPGE version of SPEX QR using Pursell method

SPEX_info SPEX_QR_PURSELL
(
    SPEX_matrix *A,            // Matrix to be factored
    SPEX_matrix **R_handle,    // upper triangular matrix
    SPEX_matrix **Q_handle     // orthogonal triangular matrix
)
{
    SPEX_info info;
    int64_t m = A->m, n = A->n;
    ASSERT( m >= n); // A should be transposed if not true
    ASSERT( A != NULL);
    // Only dense for now
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_DENSE);
    
    // Indices
    int64_t i, j, k;
    
    // A_transpose, will be overwritten with Q
    SPEX_matrix *A_T;
    // A2 = A_T A // Will be overwritten with R
    SPEX_matrix *A2;
    
    // Allocate A_T
    SPEX_CHECK(SPEX_matrix_allocate(&A_T, SPEX_DENSE, SPEX_MPZ, n, m, n*m,
        false, true, NULL));
    
    // Allocate A2
    SPEX_CHECK(SPEX_matrix_allocate(&A2, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL));
    
    // Compute A_T
    
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            // A'(i,j) = A(j,i)
            SPEX_CHECK(SPEX_mpz_set( SPEX_2D(A_T, i, j, mpz),
                          SPEX_2D(A,   j, i, mpz)));
        }
    }
    
    // Compute A2 = A'*A
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            // Compute A2(i,j) = A_T(i,:) dot A(:,j)
            for (k = 0; k < m; k++)
            {
                SPEX_CHECK(SPEX_mpz_addmul( SPEX_2D(A2,  i, j, mpz),
                                 SPEX_2D(A_T, i, k, mpz),
                                 SPEX_2D(A,   k, j, mpz)));
            }
        }
    }
    
    // Now, the factorization is computed by performing IPGE on [A2 AT]
    // We can either construct this matrix directly or do it more efficiently
    
    // Compute R and Q directly
    // A2 = A'*A will be overwriten with R (dimension n*n)
    // A' will be overwritten with Q (dimension n*m)
    
    // Perform IPGE on A2
    for (k = 0; k < n; k++)
    {
        for (i = k+1; i < n; i++)
        {
            for (j = k+1; j < n; j++)
            {
                if (i > j)
                {
                    // i > j -> entries in lower triangular portion, set them to 0
                    SPEX_CHECK(SPEX_mpz_set_ui( SPEX_2D(A2, i, j, mpz), 0));
                }
                else
                {
                    // A2(i,j) = A2(i,j)*A2(k,k)
                    SPEX_CHECK(SPEX_mpz_mul ( SPEX_2D(A2, i, j, mpz),
                                   SPEX_2D(A2, i, j, mpz),
                                   SPEX_2D(A2, k, k, mpz)));
                    // A2(i,j) = A2(i,j) - A2(i,k)*A2(i,k)
                    SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(A2, i, j, mpz),
                                     SPEX_2D(A2, k, i, mpz),
                                     SPEX_2D(A2, k, j, mpz)));
                    if (k > 0)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(A2, i, j, mpz),
                                           SPEX_2D(A2, i, j, mpz),
                                           SPEX_2D(A2, k-1, k-1, mpz)));
                    }
                }
            }
        }
        // Now that we have A2, we can compute Q
        for (i = k+1; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                SPEX_CHECK(SPEX_mpz_mul( SPEX_2D(A_T, i, j, mpz),
                              SPEX_2D(A2,  k, k, mpz),
                              SPEX_2D(A_T, i, j, mpz)));
                SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(A_T, i, j, mpz),
                                 SPEX_2D(A2,  k, i, mpz),
                                 SPEX_2D(A_T, k, j, mpz)));
                if (k > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(A_T, i, j, mpz),
                                       SPEX_2D(A_T, i, j, mpz),
                                       SPEX_2D(A2, k-1, k-1, mpz)));
                }
            }
        }
    }
    
    for (k = 1; k < n; k++)
        SPEX_CHECK(SPEX_mpz_set_ui( SPEX_2D(A2, k, 0, mpz), 0));
                                    
    (*Q_handle) = A_T;
    (*R_handle) = A2;
    return SPEX_OK;    
}
