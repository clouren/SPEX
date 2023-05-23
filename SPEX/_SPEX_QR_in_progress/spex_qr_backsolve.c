//------------------------------------------------------------------------------
// SPEX_QR/spex_qr_backsolve: Solve R x = Q' b for QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Christopher Lourenco,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_ALL ;

#include "spex_qr_internal.h"


/* Purpose: This function solves the linear system R x = Q' b for QR
 * factorization. On input, x contains the rhs. On output, x contains the exact solution
 * of the system Ax = (det A)*b R is the right triangular REF QR factor of
 * A. It is not modified on input/output
 */
 
SPEX_info spex_qr_backsolve
(
    // Output
    SPEX_matrix x,          // Solution vector to A x = det(A) * b
    // Input
    const SPEX_matrix R     // The upper triangular matrix
)
{

    SPEX_info info;
    // All inputs have been checked by the caller, asserts are
    // here as a reminder
    ASSERT(L->type == SPEX_MPZ);
    ASSERT(L->kind == SPEX_CSC);
    ASSERT(x->type == SPEX_MPZ);
    ASSERT(x->kind == SPEX_DENSE);

    int64_t k, p, j, n = R->n;
    int sgn, sgn2;

     // Set x = b
    SPEX_matrix_copy(&x, SPEX_DENSE, SPEX_MPZ, b, NULL);
    // Scale x by determinant of A'*A (R(n,n))
    for (i = 0; i < x->m*x->n; i++)
    {
        SPEX_mpz_mul( x->x.mpz[i], x->x.mpz[i], SPEX_2D(R, R->n-1, R->n-1, mpz));
    }

    // Iterate across the RHS vectors
    for (k = 0; k < x->n; k++)
    {
        // Iterate across the rows of x
        for (j = n-1; j >= 0; j--)
        {

            // Iterate across column j of L
            for (i = R->p[j]; i < R->p[j+1]-1; i++)
            {
                // If either x[i,k] or R[i,k] is 0, skip the operation
                //SPEX_MPZ_SGN(&sgn, SPEX_2D(x, R->i[i], k, mpz));
		SPEX_MPZ_SGN(&sgn, SPEX_2D(x, j, k, mpz));
                SPEX_MPZ_SGN(&sgn2, R->x.mpz[i]);
                if (sgn == 0 || sgn2 ==0 ) continue;

                // Compute x[i,k] = x[i,k] - L[i,k]*x[j,k]
                SPEX_MPZ_SUBMUL(SPEX_2D(x,  R->i[i], k, mpz),
                                R->x.mpz[i], SPEX_2D(x, j, k, mpz));

            }

            // Obtain x[j]=x[j]/R[n,n]
            SPEX_MPZ_DIVEXACT(SPEX_2D(x, j, k, mpz),
                                          SPEX_2D(x, j, k, mpz),
                                          R->x.mpz[R->p[j+1]-1]);

        }

    }
    return SPEX_OK;
}
