//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_ltsolve: Solve the system L' x = b for Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_ALLOCATION     \
    SPEX_matrix_free(&x, NULL);  \
    return SPEX_OUT_OF_MEMORY    \
    //TOASK why the error message?

#include "spex_chol_internal.h"


/* Purpose: This function solves the linear system L' x = x for Cholesky
 * factorization. On input, x contains the forward substitution solution vector
 * (that is the solution of L D x = b. On output, x contains the exact solution
 * of the system Ax = (det A)*b L is the lower triangular REF Cholesky factor of
 * A. It is not modified on input/output
 */
SPEX_info spex_Chol_backward_sub 
(
    // Output
    SPEX_matrix* x,         // Solution vector to A x = det(A) * b
    // Input
    const SPEX_matrix* L    // The lower triangular matrix
)
{
    SPEX_info info;
    // All inputs have been checked by the caller, so no input checking here
    // (otherwise it's impossible to 100% test cover)
    //No null input checking 

    // Check input
    ASSERT(L->type == SPEX_MPZ);
    ASSERT(L->kind == SPEX_CSC);
    ASSERT(x->type == SPEX_MPZ);
    ASSERT(x->kind == SPEX_DENSE);
    
    int64_t k, p, j, n = L->n;
    int sgn;
    // Iterate across the RHS vectors
    for (k = 0; k < x->n; k++)
    {
        // Iterate across the rows of x
        for (j = n-1; j >= 0; j--)
        {
            // if x[j,k] == 0 skip this operation
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, SPEX_2D(x, j, k, mpz)));
            if (sgn == 0) continue;
            
            for (p = L->p[j]+1; p < L->p[j+1]; p++)
            {
                // Compute x[j,k] = x[j,k] - L[p,k]*x[p,k]
                SPEX_CHECK( SPEX_mpz_submul( SPEX_2D(x, j, k, mpz), L->x.mpz[p], 
                                      SPEX_2D( x, L->i[p], k, mpz)));
            }
            // Compute x[j,k] = x[j,k] / L[j,j]
            SPEX_CHECK( SPEX_mpz_divexact( SPEX_2D(x, j, k, mpz), 
                        SPEX_2D(x, j, k, mpz), L->x.mpz[ L->p[j]]));
        }
    }
    return SPEX_OK;
}
