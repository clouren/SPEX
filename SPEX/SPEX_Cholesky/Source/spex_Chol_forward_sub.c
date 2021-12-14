///------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_forward_sub: Solve the system LDx = b
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF forward substitution This is
 * essentially the same as the sparse REF triangular solve applied to each
 * column of the right hand side vectors. Like the normal one, this
 * function expects that the vector x is dense. As a result,the nonzero
 * pattern is not computed and each nonzero in x is iterated across.
 * The system to solve is LDx = x
 *
 * On output, the matrix x structure is modified
 *
 */

#define SPEX_FREE_WORKSPACE      \
    SPEX_matrix_free(&h, NULL);  \

#define SPEX_FREE_ALLOCATION     \
    SPEX_FREE_WORKSPACE          \
    SPEX_matrix_free(&x, NULL);  \

#include "spex_chol_internal.h"

SPEX_info spex_Chol_forward_sub
(
    // Input/Output
    SPEX_matrix* x,              // Right hand side matrix. 
                                 // On input: contains b
                                 // On output: contains the solution of LD x = x
    // Input
    const SPEX_matrix* L,        // REF Cholesky factor of A (lower triangular)
    const SPEX_matrix* rhos      // Sequence of pivots used in factorization
)
{
    SPEX_info info;
    int64_t  i, j, p, k, n = L->n, m, mnew;
    
    ASSERT(n >= 0);
    ASSERT(L->n == x->m);
    
    int sgn;
    
    // Build the history matrix
    SPEX_matrix *h;
    SPEX_CHECK(SPEX_matrix_allocate(&h, SPEX_DENSE, SPEX_INT64, x->m, x->n,
                                    x->nzmax, false, true, NULL));

    // initialize entries of history matrix to be -1
    int xnz = x->n * x->m;
    for (i = 0; i < xnz; i++)
    {
        h->x.int64[i] = -1;
    }

    //--------------------------------------------------------------------------
    // Iterate across each RHS vector
    //--------------------------------------------------------------------------
    for (k = 0; k < x->n; k++)
    {
        //----------------------------------------------------------------------
        // Iterate accross all nonzeros in x. Assume x is dense
        //----------------------------------------------------------------------
        for (i = 0; i < n; i++)
        {
            // p is the history value of x[i,k] that is p = h[i,k]
            p = SPEX_2D(h, i, k, int64);
            // If x[i][k] = 0, can skip operations and continue to next i
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, i, k, mpz)));
            // TODO Please propagate ^ to the other codes. USE WRAPPERS DONE
            if (sgn == 0) continue;

            //------------------------------------------------------------------
            // History Update x[i,k]
            //------------------------------------------------------------------
            if (p < i-1)
            {
                // x[i,l] = x[i,k] * rhos[i-1]
                SPEX_CHECK(SPEX_mpz_mul(SPEX_2D(x, i, k, mpz),
                                    SPEX_2D(x, i, k, mpz),rhos->x.mpz[i-1]));
    
                // Only divide by the previous pivot if it is not =1
                // (the default for the first iteration)  
                if (p > -1)
                {
                    // x[i,k] = x[i,k] / rhos[p]
                    SPEX_CHECK(SPEX_mpz_divexact(SPEX_2D(x, i, k, mpz),
                                        SPEX_2D(x, i, k, mpz), rhos->x.mpz[p]));
                }
            }

            //------------------------------------------------------------------
            // IPGE updates
            //------------------------------------------------------------------
            // Access entry L[m,i]
            for (m = L->p[i]; m < L->p[i+1]; m++)
            {
                // mnew is the row index of L[m,i]
                mnew = L->i[m];
                // skip if L[m,i] is zero
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->x.mpz[m]));
                if (sgn == 0) continue;

                if (mnew > i)
                {
                    // p is the history value of x[m,k]
                    p = SPEX_2D(h, mnew, k, int64);
                    // Check if x[m,k] is zero
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, mnew, k, mpz)));
                    if (sgn == 0)
                    {
                        // In this case x[m,k] is zero,
                        // so we compute x[m,k] = 0 - l[m,i]*x[i,k]
                        SPEX_CHECK(SPEX_mpz_submul(SPEX_2D(x, mnew, k, mpz),
                                        L->x.mpz[m], SPEX_2D(x, i, k, mpz)));
                        // x[m,k] = x[m,k]/rhos[i-1] if we are not in the first
                        // iteration (in which case rhos[i-1] = 1
                        if (i > 0)
                        {
                            SPEX_CHECK(SPEX_mpz_divexact(SPEX_2D(x, mnew, k, mpz),
                                SPEX_2D(x, mnew, k, mpz), rhos->x.mpz[i-1]));
                        }
                    }
                    else
                    {
                        // In this case, x[m,k] is not equal to zero. We first
                        // check if a history update is necessary
                        if (p < i-1)
                        {
                            // x[m,k] = x[m,k] * rhos[i-1]
                            SPEX_CHECK(SPEX_mpz_mul( SPEX_2D(x, mnew, k, mpz),
                                    SPEX_2D(x, mnew, k, mpz),rhos->x.mpz[i-1]));
                            // Divide by the history pivot if we are not in the
                            // first interation (in which case rhos[p] = 1)
                            // x[m,k] = x[m,k] / rhos[p]
                            if (p > -1)
                            {
                                SPEX_CHECK(SPEX_mpz_divexact(SPEX_2D(x, mnew, k, mpz),
                                    SPEX_2D(x, mnew, k, mpz), rhos->x.mpz[p]));
                            }
                        }
                        // x[m,k] = x[m,k] * rhos[i]
                        SPEX_CHECK(SPEX_mpz_mul( SPEX_2D(x, mnew, k, mpz),
                                    SPEX_2D(x, mnew, k, mpz),rhos->x.mpz[i]));
                        // x[m,k] = x[m,k] - l[m,i] *x[i,k]
                        SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(x, mnew, k, mpz),
                                        L->x.mpz[m], SPEX_2D(x, i, k, mpz)));
                        // Divide by the previous pivot if not in iteration 0.
                        // If at iteration 0, the pivot is equal to 1 so no
                        // division is necessary
                        // x[m,k] = x[m,k] / rhos[i-1]
                        if (i > 0)
                        {
                            SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz),
                                    SPEX_2D(x, mnew, k, mpz), rhos->x.mpz[i-1]));
                        }
                    }
                    // Update the history value of x[m,k]
                    SPEX_2D(h, mnew, k, int64) = i;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Free h memory and return success
    //--------------------------------------------------------------------------
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
