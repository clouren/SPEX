//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_QR_IPGE.c: Modified Pursell Algorithm
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2022, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs REF QR via the modified Pursell Algorithm (Algorithm 1)
 * from paper.
 */

# include "spex_qr_internal.h"

/* Perform the IPGE version of SPEX QR (aka Algorithm 1 from workpage)
 */
SPEX_info SPEX_QR_IPGE
(
    SPEX_matrix A,            // Matrix to be factored
    SPEX_matrix *R_handle,    // Null on input, contains R on output
    SPEX_matrix *Q_handle     // Null on input, contains Q on output
)
{
    SPEX_info info;
    int64_t m = A->m, n = A->n;
    ASSERT( m >= n); // A should be transposed if not true
    if (m < n)
        return SPEX_PANIC;
    ASSERT( A != NULL);
    // Only dense for now
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_DENSE);

    // Indices
    int64_t i, j, k;

    // Final matrices Q and R
    SPEX_matrix Q, R;

    // Allocate R. We are performing the Thin REF QR factorization so
    // R is n*n
    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL));

    // Set Q = A
     SPEX_CHECK(SPEX_matrix_copy(&Q, SPEX_DENSE, SPEX_MPZ, A, NULL));

     // Perform Factorization
    for (k = 0; k < n; k++)
    {
        // Compute row k of R
        for (j = k; j < n; j++)
        {
            // R(k,j) = Q(:,k) dot A(:,j)
            // This is very easily parallelized
            SPEX_CHECK(SPEX_dense_mat_dot(Q, k, A, j, SPEX_2D(R,k,j,mpz)));
        }

        // IPGE update Q
        for (i = k+1; i < n; i++)
        {
            for ( j = 0; j < m; j++)
            {
                // Q(j,i) = Q(j,i)*R(k,k)
                SPEX_MPZ_MUL( SPEX_2D(Q, j, i, mpz),
                              SPEX_2D(Q, j, i, mpz),
                              SPEX_2D(R, k, k, mpz));
                // Q(j,i) = Q(j,i) - R(k,i)*Q(j,k)
                SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(Q, j, i, mpz),
                                 SPEX_2D(R, k, i, mpz),
                                 SPEX_2D(Q, j, k, mpz)));
                if (k > 0)
                {
                    // Q(j,i) = Q(j,i)/R(k-1,k-1)
                    SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(Q, j, i, mpz),
                                       SPEX_2D(Q, j, i, mpz),
                                       SPEX_2D(R, k-1, k-1, mpz)));
                }
            }
        }
    }

    (*Q_handle) = Q;
    (*R_handle) = R;
    return SPEX_OK;
}
