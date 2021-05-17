///------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_forward_sub: Solve the system LDx = b
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

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

#define SPEX_FREE_WORKSPACE            \
{                                      \
SPEX_matrix_free(&h, NULL);             \
}

#include "spex_chol_internal.h"

SPEX_info spex_Chol_forward_sub
(
    const SPEX_matrix *L,   // lower triangular matrix
    SPEX_matrix *x,         // right hand side matrix of size n*numRHS
    const SPEX_matrix *rhos // sequence of pivots used in factorization
)
{
    SPEX_info info;
    int64_t  i, j, p, k, n, m, mnew;
    // Size of x vector
    n = L->n;

    ASSERT(n >=0)
    
    int sgn;
    // calloc is used, so that h is initialized for SPEX_FREE_WORKSPACE
    
    // Build the history matrix
    SPEX_matrix *h;
    SPEX_CHECK (SPEX_matrix_allocate(&h, SPEX_DENSE, SPEX_INT64, x->m, x->n,
        x->nzmax, false, true, NULL));

    // initialize entries of history matrix to be -1
    for (i = 0; i < x->nzmax; i++)
    {
        h->x.int64[i] = -1;
    }
            
    //--------------------------------------------------------------------------
    // Iterate across each RHS vector
    //--------------------------------------------------------------------------

    for (k = 0; k < x->n; k++)
    {
        //----------------------------------------------------------------------
        // Iterate accross all nonzeros int64_t* x. Assume x is dense
        //----------------------------------------------------------------------
        for (i = 0; i < n; i++)
        {
            p = SPEX_2D(h, i, k, int64);
            // If x[i][k] = 0, can skip operations and continue to next i
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, i, k, mpz)));
            if (sgn == 0) {continue;}

            //------------------------------------------------------------------
            // History Update
            //------------------------------------------------------------------
            if (p < i-1)
            {
                // x[i] = x[i] * rhos[i-1]
                SPEX_CHECK(SPEX_mpz_mul( SPEX_2D(x, i, k, mpz), SPEX_2D(x, i, k, mpz),
                                 rhos->x.mpz[i-1]));
                // x[i] = x[i] / rhos[p]
                if (p > -1)
                {
                    SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(x, i, k, mpz), SPEX_2D(x, i, k, mpz),
                                          rhos->x.mpz[p]));
                }
            }

            //------------------------------------------------------------------
            // IPGE updates
            //------------------------------------------------------------------
            // Access the Lmi
            for (m = L->p[i]; m < L->p[i+1]; m++)
            {
                // Location of Lmi
                mnew = L->i[m];
                // skip if Lx[m] is zero
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->x.mpz[m]));
                if (sgn == 0) {continue;}
                // m > i
                if (mnew > i)
                {
                    p = SPEX_2D(h, mnew, k, int64);
                    // x[mnew] is zero
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, mnew, k, mpz)));
                    if (sgn == 0)
                    {
                        // x[m] = x[m] - lmi xi
                        SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(x, mnew, k, mpz), L->x.mpz[m],
                                            SPEX_2D(x, i, k, mpz)));
                        // x[m] = x[m] / rhos[i-1]
                        if (i > 0)
                        {
                            SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz), rhos->x.mpz[i-1]));
                        }
                    }
                    else
                    {
                        // History update if necessary
                        if (p < i-1)
                        {
                            // x[m] = x[m] * rhos[i-1]
                            SPEX_CHECK(SPEX_mpz_mul( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz),
                                             rhos->x.mpz[i-1]));
                            // x[m] = x[m] / rhos[p]
                            if (p > -1)
                            {
                                SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz),
                                                      rhos->x.mpz[p]));
                            }
                        }
                        // x[m] = x[m] * rhos[i]
                        SPEX_CHECK(SPEX_mpz_mul( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz),
                                         rhos->x.mpz[i]));
                        // x[m] = x[m] - lmi xi
                        SPEX_CHECK(SPEX_mpz_submul( SPEX_2D(x, mnew, k, mpz), L->x.mpz[m], SPEX_2D(x, i, k, mpz)));
                        // x[m] = x[m] / rhos[i-1]
                        if (i > 0)
                        {
                            SPEX_CHECK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz), 
                                                  rhos->x.mpz[i-1]));
                        }
                    }
                    SPEX_2D(h, mnew, k, int64) = i;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Free h memory
    //--------------------------------------------------------------------------
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}

