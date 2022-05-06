//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Solve: Solves the system LDL' x = b.
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE        \
{                                  \
    SPEX_matrix_free(&b2, option); \
    SPEX_matrix_free(&x2, option);  \
}

# define SPEX_FREE_ALL             \
{                                  \
    SPEX_FREE_WORKSPACE            \
    SPEX_matrix_free(&x, NULL);    \
}
//TODO have a separate "clusterVersion" for all the prints, so they are here
#include "spex_chol_internal.h"
    
/* Purpose: This function solves the linear system LDL' x = b.
 *
 * Input arguments:
 * 
 * x_handle:        A handle to the solution matrix. NULL on input.
 *                  On output contains a pointer to the solution vector(s).
 * 
 * F:               Factorization struct. Contains the REF Cholesky factorization
 *                  of A along with permutation and elimination tree information.
 * 
 * b:               Right hand side vector(s). Must be SPEX_MPZ and SPEX_DENSE.
 * 
 * option:          Command options. Default if NULL.
 */

SPEX_info SPEX_Chol_solve       // Solves the linear system LDL' x = b
(
    // Output
    SPEX_matrix** x_handle,           // On input: NULL.
                                      // On output: Exact rational solution (SPEX_MPQ)
                                      // to the system. 
    // Input
    const SPEX_factorization* F,      // Cholesky factorization. Must contain the
                                      // elimination tree of A,  
                                      // column pointers of L, exact number of
                                      // nonzeros of L and permutation used //TOCHECK
    const SPEX_matrix* b,             // Right hand side vector(s).
                                      // Must be SPEX_MPZ and SPEX_DENSE.
    const SPEX_options* option        // Command options. Default if NULL.
)
{
    SPEX_info info;

    if (!spex_initialized()) return SPEX_PANIC; //TODO use this one in the other functions
  
    // Check the inputs
    ASSERT(F->kind==SPEX_CHOLESKY_FACTORIZATION); //TODO make part of if!

    if (!x_handle || b->type != SPEX_MPZ || b->kind != SPEX_DENSE)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
  
    // Indices, m, n and nz
    int64_t i, j, m = b->m, n = b->n, nz;
    // Set the number of nonzeros in b
    nz = m * n;
  
    // det is the determinant of the PAP matrix. It is obtained for free
    // from the SPEX Cholesky factorization det = rhos[L->n-1] = L[L->n,L->n]
    mpz_t* det = NULL;
    det = &(F->rhos->x.mpz[F->L->n-1]);
    size_t bitsDet = mpz_sizeinbase((*det),2); //CLUSTER
    //gmp_printf("%llu, ", (bitsDet)); //CLUSTER

    // Declare workspace and output
    // x2 is the (permuted) final solution vector
    SPEX_matrix *x2 = NULL;
    // x is the final solution vector returned to user
    SPEX_matrix *x = NULL;
    // b2 is the permuted right hand side vector(s)
    SPEX_matrix *b2 = NULL;

    // Allocate memory for b2
    SPEX_CHECK( SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPZ, m, n, 
                                     nz, false, true, option) );

    // Set b2[i,j] = b[pinv[i],j] to permute the RHS vectors. We
    // will then solve PAP b2 = b2 overwriting these permuted vectors
    // with the scaled integer solution
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)   
        {
            SPEX_CHECK( SPEX_mpz_set(SPEX_2D(b2, F->Pinv_perm[i], j, mpz),
                                     SPEX_2D(b, i, j, mpz)) );
        }
    }

    // Forward substitution, b2 = L \ b2. Note that b2 is overwritten    
    SPEX_CHECK( spex_chol_forward_sub(b2, F->L, F->rhos) );

    // Apply the determinant to b2: b2 = b2*det(A)
    for (i = 0; i < nz; i++)
    {
        SPEX_CHECK( SPEX_mpz_mul(b2->x.mpz[i], b2->x.mpz[i], (*det)) );
    }
    
    // Backsolve, b2 = L' \ b2. Note that, again, b2 is overwritten
    SPEX_CHECK( spex_chol_backward_sub(b2, F->L) );

    // Allocate memory for the (permuted) solution vector x2 = b2/det
    SPEX_CHECK( SPEX_matrix_allocate(&x2, SPEX_DENSE, SPEX_MPQ, m, n,
                                     0, false, true, option) ) ;
/* TODO DELETE ME
    size_t bitsSmall=mpz_sizeinbase(b2->x.mpz[0],2); //CLUSTER
    size_t bitsLarge=0; //CLUSTER
    size_t bitsTemp; //CLUSTER
    int64_t bitsSum=0; //CLUSTER
    // Set x2[i] = b2[i]/det
*/
    // Compute the solution vector to the system PAP x2 = b2
    // This is done by setting x2 = b2/det
    for (i = 0; i < nz; i++)
    {
        /* TODO DELETE ME
        bitsTemp=mpz_sizeinbase(b2->x.mpz[i],2); //CLUSTER
        if(bitsTemp>bitsLarge) //CLUSTER
        { 
            bitsLarge=bitsTemp; //CLUSTER
        }
        else if(bitsTemp<bitsSmall) //CLUSTER
        {
                bitsSmall=bitsTemp; //CLUSTER
        }
        bitsSum=bitsSum+(int64_t)bitsTemp; //CLUSTER
*/
        // Set the numerator of x2[i] = b2[i]
        SPEX_CHECK( SPEX_mpq_set_num(x2->x.mpq[i], b2->x.mpz[i]) );
        // Set the denominator of x2[i] = det
        SPEX_CHECK( SPEX_mpq_set_den(x2->x.mpq[i], (*det)) );
        // Remove common factors from numerator and denominator
        SPEX_CHECK( SPEX_mpq_canonicalize(x2->x.mpq[i]) );
    }
    /* TODO DELETE ME
    gmp_printf("%llu, %llu, %llu, ", bitsSmall, bitsLarge, bitsSum); //CLUSTER
    */
    // Allocate memory for x which is the final solution vector of the 
    // original (unpermuted) linear system. 
    SPEX_CHECK( SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPQ, m, n, 0, false, true, option) );
    
    // Undo the inverse row permutation on x2. This is done by applying the row permutation
    // to the indices of x. That is x[p[i], j] = x2[i,j]
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            SPEX_CHECK( SPEX_mpq_set(SPEX_2D(x, F->P_perm[i], j, mpq), 
                                     SPEX_2D(x2, i, j, mpq)) );
        }
    }

    //--------------------------------------------------------------------------
    // Scale the solution if necessary.
    //--------------------------------------------------------------------------
    // When the forward/backsolve is complete, the entries in
    // x are rational, but are solving the scaled linear system
    // A' x = b' (that is if A had input which was rational or floating point
    // and had to be converted to integers). Thus, the scale here is used
    // to convert x into into the actual solution of A x = b. 
    SPEX_CHECK( SPEX_scale(x, F->scale_for_A, b->scale, option) );

    // Set output, free memory, return success
    (*x_handle) = x;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
