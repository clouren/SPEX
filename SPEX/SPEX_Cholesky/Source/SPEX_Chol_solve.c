//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Solve: Solve the SPD linear system after factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE        \
{                                  \
    SPEX_matrix_free(&b2, option); \
}

# define SPEX_FREE_ALL             \
{                                  \
    SPEX_FREE_WORKSPACE            \
    SPEX_matrix_free(&x, NULL);    \
}
//TODO maybe have a separate "clusterVersion" for all the prints, so they are def not in the release version
#include "spex_chol_internal.h"

/* Purpose: This function solves the linear system LDL' x = b.
 *
 * Input arguments:
 * 
 * x_handle:        A handle to the solution matrix. On input this is NULL,
 *                  on output x_handle contains a pointer to the solution
 *                  vector(s)
 * 
 * PAP:             Permuted version of the input matrix
 * 
 * A:               Nonpermuted (original) input matrix
 * 
 * b:               Right hand side vector(s)
 * 
 * rhos:            Sequence of factorization pivots
 * 
 * L:               Lower triangular matrix
 * 
 * S:               Symbolic analysis struct (contains row/column permutation
 *                  and its inverse
 * 
 * option:          Command options
 * 
 */

SPEX_info SPEX_Chol_solve
(
    // Output
    SPEX_matrix** x_handle,           // On input: undefined.
                                      // On output: Rational solution (SPEX_MPQ)
                                      // to the system. 
    // input/output:
    SPEX_factorization *F,  // The non-updatable Cholesky factorization.
                            // Mathematically, F is unchanged.  However, if F
                            // is updatable on input, it is converted to
                            // non-updatable.  If F is already non-updatable,
                            // it is not modified.
    // input:
    const SPEX_matrix* b,             // Right hand side vector
    const SPEX_options* option        // command options
)
{
    SPEX_info info;

    if (!spex_initialized()) return SPEX_PANIC;
  
    // Check the inputs
    ASSERT(!x_handle);
    ASSERT(b->type == SPEX_MPZ);
    ASSERT(b->kind == SPEX_DENSE);
    ASSERT(F->kind==SPEX_CHOLESKY_FACTORIZATION);

    if (!x_handle || b->type != SPEX_MPZ || b->kind != SPEX_DENSE)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    // convert the factorization F to non-updatable
    info = SPEX_factorization_convert(F, false, option);
    if (info != SPEX_OK) return info;
    
    //int64_t i, j, nz, n = F->L->n;

    // det is the determinant of the PAP matrix. It is obtained for free
    // from the SPEX Cholesky factorization det = rhos[n-1] = L[n,n]
    mpz_t* det = NULL;

    //--------------------------------------------------------------------------
    // Declare workspace and output
    //--------------------------------------------------------------------------
    // x is the permuted final solution vector returned to the user
    SPEX_matrix *x = NULL;
    // b2 is the permuted right hand side vector(s)
    SPEX_matrix *b2 = NULL;

    //--------------------------------------------------------------------------
    // get b2 = Pinv*b
    //--------------------------------------------------------------------------

    SPEX_CHECK (spex_permute_dense_matrix (&b2, b, F->Pinv_perm, option)) ;

    //--------------------------------------------------------------------------
    // Forward substitution, b2 = L \ b2. Note that b2 is overwritten    
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_chol_forward_sub(b2, F->L, F->rhos));

    //--------------------------------------------------------------------------
    // Apply the determinant to b2, b2 = det*b2
    //--------------------------------------------------------------------------

    // Set the value of the determinant det = rhos[n-1] 
    det = &(F->rhos->x.mpz[F->L->n-1]);
    // FIXME unused variable (commented by Jinhao)
    // size_t bitsDet = mpz_sizeinbase((*det),2); //CLUSTER
    //gmp_printf("%llu, ", (bitsDet)); //CLUSTER
    
    SPEX_CHECK(spex_matrix_mul(b2, (*det) ));
    
    //--------------------------------------------------------------------------
    // Backsolve, b2 = L' \ b2. Note that, again, b2 is overwritten
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_chol_backward_sub(b2, F->L));

    //--------------------------------------------------------------------------
    // get real solution x by applying both permutation and scale
    // x = P*b2/scale
    //--------------------------------------------------------------------------    // Scale is the scaling factor for the solution vectors.
    // When the forward/backsolve is complete, the entries in
    // x/det are rational, but are solving the scaled linear system
    // A' x = b' (that is if A had input which was rational or floating point
    // and had to be converted to integers). Thus, the scale here is used
    // to convert x into into the actual solution of A x = b. 
    // Mathematically, set scale = b->scale * rhos[n-1] / PAP->scale
    SPEX_CHECK(SPEX_mpq_set_z(b2->scale, (*det)));
    SPEX_CHECK(SPEX_mpq_mul(b2->scale, b2->scale, b->scale));
    SPEX_CHECK(SPEX_mpq_div(b2->scale, b2->scale, F->scale_for_A));

    // allocate space for x as dense MPQ matrix
    SPEX_CHECK (SPEX_matrix_allocate (&x, SPEX_DENSE, SPEX_MPQ, b->m, b->n,
        0, false, true, option));
/*
    size_t bitsSmall=mpz_sizeinbase(b2->x.mpz[0],2); //CLUSTER
    size_t bitsLarge=0; //CLUSTER
    size_t bitsTemp; //CLUSTER
    int64_t bitsSum=0; //CLUSTER
    // Set x[i] = b2[i]/det
*/

    // obtain x from permuted b2 with scale applied
    for (int64_t i = 0 ; i < b->m ; i++)
    {
        int64_t pi = F->P_perm[i];
        for (int64_t j = 0 ; j < b->n ; j++)
        {
            /*
            bitsTemp=mpz_sizeinbase(SPEX_2D(b2,  i, j, mpz),2); //CLUSTER
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
            SPEX_CHECK(SPEX_mpq_set_z(SPEX_2D(x,  pi, j, mpq),
                                      SPEX_2D(b2,  i, j, mpz)));
            SPEX_CHECK(SPEX_mpq_div(SPEX_2D(x,  pi, j, mpq),
                                    SPEX_2D(x,  pi, j, mpq), b2->scale));
        }
    }

    /*
    gmp_printf("%llu, %llu, %llu, ", bitsSmall, bitsLarge, bitsSum); //CLUSTER
    */

    // Set output, free memory
    (*x_handle) = x;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
