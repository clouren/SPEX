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
    for (i = 0; i < n; i++)            \
    {                                  \
        SPEX_FREE(h[i]);               \
    }                                  \
    SPEX_FREE(h);                      \
}

#include "spex_chol_internal.h"

SPEX_info spex_Chol_forward_sub
(
    SPEX_matrix *L,         // lower triangular matrix
    SPEX_matrix *x,         // right hand side matrix of size n*numRHS
    SPEX_matrix *rhos       // sequence of pivots used in factorization
)
{
    SPEX_info ok;
    int64_t  i, j, p, k, n, m, mnew, **h;
    // Size of x vector
    n = L->n;

    int sgn;
    // calloc is used, so that h is initialized for SPEX_FREE_WORKSPACE
    h = (int64_t**) SPEX_calloc(n, sizeof(int64_t*));
    if (!h)
    {
        return SPEX_OUT_OF_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        h[i] = (int64_t*) SPEX_malloc(x->n* sizeof(int64_t));
        if (!h[i])
        {
            SPEX_FREE_WORKSPACE;
            return SPEX_OUT_OF_MEMORY;
        }
        for (j = 0; j < x->n; j++)
        {
            h[i][j] = -1;
        }
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
            p = h[i][k];
            // If x[i][k] = 0, can skip operations and continue to next i
            OK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, i, k, mpz)));
            if (sgn == 0) {continue;}

            //------------------------------------------------------------------
            // History Update
            //------------------------------------------------------------------
            if (p < i-1)
            {
                // x[i] = x[i] * rhos[i-1]
                OK(SPEX_mpz_mul( SPEX_2D(x, i, k, mpz), SPEX_2D(x, i, k, mpz),
                                 rhos->x.mpz[i-1]));
                // x[i] = x[i] / rhos[p]
                if (p > -1)
                {
                    OK(SPEX_mpz_divexact( SPEX_2D(x, i, k, mpz), SPEX_2D(x, i, k, mpz),
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
                OK(SPEX_mpz_sgn(&sgn, L->x.mpz[m]));
                if (sgn == 0) {continue;}
                // m > i
                if (mnew > i)
                {
                    p = h[mnew][k];
                    // x[mnew] is zero
                    OK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, mnew, k, mpz)));
                    if (sgn == 0)
                    {
                        // x[m] = x[m] - lmi xi
                        OK(SPEX_mpz_submul( SPEX_2D(x, mnew, k, mpz), L->x.mpz[m],
                                            SPEX_2D(x, i, k, mpz)));
                        // x[m] = x[m] / rhos[i-1]
                        if (i > 0)
                        {
                            OK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz), rhos->x.mpz[i-1]));
                        }
                    }
                    else
                    {
                        // History update if necessary
                        if (p < i-1)
                        {
                            // x[m] = x[m] * rhos[i-1]
                            OK(SPEX_mpz_mul( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz),
                                             rhos->x.mpz[i-1]));
                            // x[m] = x[m] / rhos[p]
                            if (p > -1)
                            {
                                OK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz),
                                                      rhos->x.mpz[p]));
                            }
                        }
                        // x[m] = x[m] * rhos[i]
                        OK(SPEX_mpz_mul( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz),
                                         rhos->x.mpz[i]));
                        // x[m] = x[m] - lmi xi
                        OK(SPEX_mpz_submul( SPEX_2D(x, mnew, k, mpz), L->x.mpz[m], SPEX_2D(x, i, k, mpz)));
                        // x[m] = x[m] / rhos[i-1]
                        if (i > 0)
                        {
                            OK(SPEX_mpz_divexact( SPEX_2D(x, mnew, k, mpz), SPEX_2D(x, mnew, k, mpz), 
                                                  rhos->x.mpz[i-1]));
                        }
                    }
                    h[mnew][k] = i;
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

