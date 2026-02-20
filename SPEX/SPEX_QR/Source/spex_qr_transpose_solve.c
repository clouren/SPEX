//------------------------------------------------------------------------------
// SPEX_QR/Source/spex_qr_transpose_solve.c: Solve exactly x = Q D * (R^T D \ b)
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2026, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function solves the sparse x = Q D ( R^T D \ b).
 * Q D R are the REF QR factorization of A^T and thus this function is intended
 * for use when solving Ax = b and A is rectangular with more columns than rows.
 *
 * This function first solves y = R^T D \b and then calculates x as
 * x = Q D y
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Unitialized on input.
 *           on output, contains the exact rational solution of the system
 *
 * b:        Set of RHS vectors
 *
 * F:        QR factorization of A^T.
 *
 * option:   command options
 */

#define SPEX_FREE_WORKSPACE               \
    {                                     \
        SPEX_matrix_free(&b_new, option); \
        SPEX_free(Qinv_perm);             \
    }

#define SPEX_FREE_ALL                 \
    {                                 \
        SPEX_FREE_WORKSPACE           \
        SPEX_matrix_free(&x, option); \
    }

#include "spex_qr_internal.h"
#include "spex_lu_internal.h"

SPEX_info spex_qr_transpose_solve(
    // Output
    SPEX_matrix *x_handle, // On input: undefined.
                           // On output: Rational solution (SPEX_MPQ)
                           // to the system.
    // input
    const SPEX_factorization F, // The QR factorization.
    const SPEX_matrix b,        // Right hand side vector
    const SPEX_options option   // command options
)
{
    SPEX_info info;
    // Check inputs, the number of columns of Q must equal the number of rows of b
    // Also, both matrices must be mpz and dense
    // Ensure SPEX is initialized
    if (!spex_initialized())
    {
        return SPEX_PANIC;
    }

    // Check the inputs
    if (!x_handle || b->type != SPEX_MPZ || b->kind != SPEX_DENSE || F->kind != SPEX_QR_FACTORIZATION)
    {
        return SPEX_INCORRECT_INPUT;
    }

    // Step 1 is a transpose triangular solve R^T D y = b
    // A can be rank deficient so if A is rank deficient R contains n-rank rows of zeros
    // meanign that R^T contains n-rank columns of zeros. So the transpose solve should just loop
    // through columns 0 to n-rank and do the typical forward sub.


    // Declare x and b_new
    SPEX_matrix x = NULL, b_new = NULL, b2 = NULL;

    int64_t *Qinv_perm = NULL;
    int64_t i, j, p, k;

    // Permute b_new

    // TODO Check this. We have AT P = Q D R, so A = P R^T D Q^T
    // Thus, P R^T D Q^T x = b so we have bnew = P^T b. Since here
    // P is F->Q_perm, I think to match LU/Cholesky here F->Q_perm goes

    SPEX_CHECK (spex_permute_dense_matrix (&b_new, b, F->Q_perm, option));

    // Transpose forward solve. Set b_new = (R^T D) \ b_new
    SPEX_CHECK( spex_qr_transpose_forward_sub( F->R, F->rank, b_new, F->rhos));

    // Now we have b_new = (R' D) \ b. The next step is to calculate x
    // x = Q D b_new
    // We will do a scaling with D first and then do the dot products with Q

    // Loop through each RHS vector
    for (k = 0; k < b_new->n; k++)
    {
        // Need to multiply each entry by the associated entry in D
        // Recall that D[j,j] = rhos[j]*rhos[j-1]
        // In order to avoid an if in the inner loop, we will do b[0]
        // here and then the rest in the for
        SPEX_CHECK( SPEX_mpz_mul( SPEX_2D(b_new, 0, k, mpz),
                                  SPEX_2D(b_new, 0, k, mpz), F->rhos->x.mpz[0]));
        // Only the entries in b_new[0..rank] are nonzero. Loop through
        // what's left
        for (j = 1; j < F->rank; j++)
        {
            // Compute D[j,j] * b_new[j]
            // First b_new[j] = b_new[j]*rhos[j-1]
            SPEX_CHECK( SPEX_mpz_mul( SPEX_2D(b_new, j, k, mpz),
                                      SPEX_2D(b_new, j, k, mpz),
                                      F->rhos->x.mpz[j-1]));

            SPEX_CHECK( SPEX_mpz_mul( SPEX_2D(b_new, j, k, mpz),
                                      SPEX_2D(b_new, j, k, mpz),
                                      F->rhos->x.mpz[j]));
        }
    }

    // Now, b_new = D*b_new
    // All that's left is to calculate Q*b_new
    // Now, we have to compute Q*(D b_new).

    // We need b2 for the final multiplication
    SPEX_CHECK(SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPZ, b->m, b->n, 0,
                                    false, true, NULL));

    // Loop through each RHS vector
    for (k = 0; k < b_new->n; k++)
    {
        // Loop through columns 1:rank of Q
        // If A is rank deficient Q will contain
        // n-rank columns of zeros
        for (j = 0; j < F->rank; j++)
        {
            // Loop through the nonzeros in each column
            // b2[i] += Q[i,j]*b_new[i]
            for (p = F->Q->p[j]; p < F->Q->p[j+1]; p++)
            {
                i = F->Q->i[p];
                SPEX_CHECK( SPEX_mpz_addmul( SPEX_2D(b2, i, k, mpz), SPEX_2D(b_new, i, k, mpz), F->Q->x.mpz[p]));
            }
        }
    }

    int64_t n = F->Q->n;
    int sgn;

    Qinv_perm = (int64_t *)SPEX_malloc(n * sizeof(int64_t));
    if (!Qinv_perm)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    for (k = 0; k < n; k++)
    {
        int64_t index = F->Q_perm[k];
        Qinv_perm[index] = k;
    }

    //--------------------------------------------------------------------------
    // x = b2/scale
    //--------------------------------------------------------------------------
    // set scale = b->scale / A_scale
    SPEX_MPQ_SET_Z(b2->scale, b->scale);
    SPEX_MPQ_DIV(b2->scale, b2->scale, F->scale_for_A);

    // allocate space for x as dense MPQ matrix
    SPEX_CHECK(SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPQ, F->Q->n, b->n,
                                    0, false, true, option));

    // obtain x from permuted b2 with scale applied
    for (i = 0; i < F->Q->n; i++)
    {
        int64_t qi = Qinv_perm[i];
        for (j = 0; j < b->n; j++)
        {
            SPEX_MPQ_SET_Z(SPEX_2D(x, qi, j, mpq),
                           SPEX_2D(b2, i, j, mpz));
            SPEX_MPQ_DIV(SPEX_2D(x, qi, j, mpq),
                         SPEX_2D(x, qi, j, mpq), b2->scale);
        }
    }

    //--------------------------------------------------------------------------
    // Return result and free workspace
    //--------------------------------------------------------------------------
    (*x_handle) = x;

    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
