//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_QR_PURSELL2.c: Pursell Algorithm v2
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs REF QR via the Pursell Algorithm on A'*A. This version
 * explicitly constructs [A'A A']. 
 */


# define SPEX_FREE_ALL              \
    SPEX_matrix_free(&A_T, NULL);   \
    SPEX_matrix_free(&A2, NULL);    \
    SPEX_matrix_free(&A3, NULL);    \


# include "spex_qr_internal.h"

    
/* Perform the IPGE version of SPEX QR using Pursell method
 */
SPEX_info SPEX_QR_PURSELL2
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
    
    // Final matrices Q and R
    SPEX_matrix *Q, *R;
    
    // A_transpose
    SPEX_matrix *A_T;
    // A2 = A_T A
    SPEX_matrix *A2;
    SPEX_matrix *A3;
    // Allocate R
    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL));
    
    // Allocate Q
    SPEX_CHECK(SPEX_matrix_allocate(&Q, SPEX_DENSE, SPEX_MPZ, m, n, m*n,
        false, true, NULL));
     
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
            SPEX_CHECK(SPEX_mpz_set( SPEX_2D(A_T, i, j, mpz),
                          SPEX_2D(A, j, i, mpz)));
        }
    }
    
    // Compute A2 = A_T*A
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            // Compute A2(i,j) = A_T(i,:) dot A(:,j)
            for (k = 0; k < m; k++)
            {
                SPEX_CHECK(SPEX_mpz_addmul( SPEX_2D(A2, i, j, mpz),
                                 SPEX_2D(A_T, i, k, mpz),
                                 SPEX_2D(A, k, j, mpz)));
            }
        }
    }
    
    // Now, the factorization is computed by performing IPGE on [A2 AT]
    // We can either construct this matrix directly or do it more efficiently
    
    // Allocate A3
    SPEX_CHECK(SPEX_matrix_allocate(&A3, SPEX_DENSE, SPEX_MPZ, n, n+m, n*(n+m),
        false, true, NULL));
    // Populate first n*n portion
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            SPEX_CHECK(SPEX_mpz_set( SPEX_2D(A3, i, j, mpz),
                          SPEX_2D(A2, i, j, mpz)));
        }
    }
    
    // Populate last n*m portion
    for(i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            SPEX_CHECK(SPEX_mpz_set( SPEX_2D(A3, i, j+n, mpz),
                          SPEX_2D(A_T, i, j, mpz)));
        }
    }
    
    // Now, perform IPGE on A3
    for (k = 0; k < n; k++)
    {
        for (i = k+1; i < n; i++)
        {
            //for (j = k+1; j < A3->n; j++)
            for (j = i; j < A3->n; j++)
            {
                //TODO This is slower because it's not operating on only the upper
                // triangular portion of R unlike the above code. Fix this!
                SPEX_CHECK(SPEX_mpz_mul ( SPEX_2D(A3, i, j, mpz),
                               SPEX_2D(A3, i, j, mpz),
                               SPEX_2D(A3, k, k, mpz)));
                SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(A3, i, j, mpz),
                                 SPEX_2D(A3, k, i, mpz),
                                 SPEX_2D(A3, k, j, mpz)));
                if (k > 0)
                    SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(A3, i, j, mpz),
                                       SPEX_2D(A3, i, j, mpz),
                                       SPEX_2D(A3, k-1, k-1, mpz)));
            }
        }
    }
    
    // Now get Q and R
    // Populate first n*n portion
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i <= j)
            {
            SPEX_CHECK(SPEX_mpz_set( SPEX_2D(R, i, j, mpz),
                          SPEX_2D(A3, i, j, mpz)));
            }
            else
                SPEX_CHECK(SPEX_mpz_set_ui( SPEX_2D(R,i,j,mpz), 0));
        }
    }
    
    // Populate last n*m portion
    for(i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            SPEX_CHECK(SPEX_mpz_set( SPEX_2D(Q, j, i, mpz),
                          SPEX_2D(A3, i, j+n, mpz)));
        }
    }
    
    (*Q_handle) = Q;
    (*R_handle) = R;
    SPEX_FREE_ALL;
    return SPEX_OK;
}
