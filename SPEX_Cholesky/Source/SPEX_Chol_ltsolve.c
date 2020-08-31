//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_ltsolve: Solve the system L' x = b for Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.
//------------------------------------------------------------------------------

#define FREE_WORKSPACE          \
return SPEX_OUT_OF_MEMORY;      \

#include "SPEX_Chol.h"

/* Purpose: This solves the system L'x = b for Cholesky factorization 
 * On input, L contains the lower triangular matrix. x has the solution
 * to the linear system from forward substitution
 */
SPEX_info SPEX_Chol_ltsolve 
(
    SPEX_matrix *L,     // The lower triangular matrix
    SPEX_matrix *x      // Solution vector
)
{
    SPEX_info ok;
    // Check input
    SPEX_REQUIRE(L, SPEX_CSC, SPEX_MPZ);
    SPEX_REQUIRE(x, SPEX_DENSE, SPEX_MPZ);
    
    int64_t p, j, n, k;
    // Set n
    n = L->n;
    // Iterate across the RHS vectors
    for (int64_t k = 0; k < x->n; k++)
    {
        // Iterate across the rows of x
        for (j = n-1; j >= 0; j--)
        {
            // if x(i,k) == 0 skip this operation
            if ( mpz_sgn ( SPEX_2D(x, j,k, mpz)) == 0) continue;
            
            for (p = L->p[j]+1; p < L->p[j+1]; p++)
            {
                OK ( SPEX_mpz_submul( SPEX_2D(x, j, k, mpz), L->x.mpz[p], 
                                      SPEX_2D( x, L->i[p], k, mpz)));
            }
            OK( SPEX_mpz_divexact( SPEX_2D(x, j, k, mpz), SPEX_2D(x, j, k, mpz), L->x.mpz[ L->p[j]]));
        }
    }
    return SPEX_OK;
}
