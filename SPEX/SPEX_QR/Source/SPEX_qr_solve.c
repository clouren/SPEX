//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_qr_solve.c: Solve exactly x = R \ Q^T b
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs the sparse R \ Q^T b. First it computes R(n,n)*Q^T b has already
 * been computed. Then does
 * back solve.  Returns x as mpq. //TOASK what to return x as, rn I'm mimic-ing LU where it is returned as a rational number
 */

#define SPEX_FREE_WORKSPACE        \
{                                  \
    SPEX_matrix_free(&b_perm, option); \
}

# define SPEX_FREE_ALL             \
{                                  \
    SPEX_FREE_WORKSPACE            \
    SPEX_matrix_free(&b_new, option); \
}

# include "spex_qr_internal.h"
# include "spex_lu_internal.h"

SPEX_info SPEX_qr_solve
(
    SPEX_matrix *x_handle, // Solution
    SPEX_factorization F,
    const SPEX_matrix b,        // Q^T * b
    const SPEX_options option
)
{
    SPEX_info info;
    // Check inputs, the number of columns of Q' must equal the number of rows of b
    // Also, both matrices must be mpz and dense
    ASSERT( Q->m == b->m);
    ASSERT( Q->type == SPEX_MPZ);
    ASSERT( b->type == SPEX_MPZ);
    ASSERT( Q->kind == SPEX_CSC);
    ASSERT( b->kind == SPEX_DENSE);

    SPEX_matrix b_new = NULL, b_perm=NULL, x=NULL;
    int64_t k, p, i,j;
    int64_t qi;

    // b->new has Q->n rows and b->n columns
    SPEX_CHECK(SPEX_matrix_allocate(&b_new, SPEX_DENSE, SPEX_MPZ, b->m, b->n, 0,
        false, true, NULL));
    // Need to compute b_new[i] = R(n,n)* Q'[i,:] dot b[i]
    // This is equivalent to b_new[i] = R(n,n)* Q[:,i] dot b[i]
    
    SPEX_CHECK (spex_permute_dense_matrix (&b_perm, b, F->Q_perm, option)); 
    // Iterate across every RHS vector
    for (k = 0; k < b->n; k++) //if b is a vector this will only be once, but b can be a matrix or rhs's
    {
        // Compute b[j,k]
        for(j=0;j<F->Q->n;j++)
        {
            for(p=F->Q->p[j]; p < F->Q->p[j+1]; p++)
            {
                i=F->Q->i[p];
                SPEX_MPZ_ADDMUL(SPEX_2D(b_new, j, k, mpz),F->Q->x.mpz[p],SPEX_2D(b_perm, i, k, mpz));
            }
            //F->rhos->x.mpz[F->R->n-1] is the determinant
            SPEX_MPZ_MUL (SPEX_2D(b_new, j, k, mpz),SPEX_2D(b_new, j, k, mpz),F->rhos->x.mpz[F->R->n-1]);
        }
    }
    
    //backwards substitution
    //Solves Rx=b_new (overwrites b_new into x)
    SPEX_CHECK (spex_left_lu_back_sub(F->R,b_new));
    
    //--------------------------------------------------------------------------
    // x = Q*b_new/scale
    //--------------------------------------------------------------------------
    // set scale = b->scale * rhos[n-1] / A_scale
    SPEX_MPQ_SET_Z(b_new->scale, F->rhos->x.mpz[F->R->n-1]);
    SPEX_MPQ_MUL(b_new->scale, b_new->scale, b->scale);
    SPEX_MPQ_DIV(b_new->scale, b_new->scale, F->scale_for_A);

    // allocate space for x as dense MPQ matrix
    SPEX_CHECK (SPEX_matrix_allocate (&x, SPEX_DENSE, SPEX_MPQ, b->m, b->n,
        0, false, true, option));
    
    // obtain x from permuted b2 with scale applied
    for (i = 0 ; i < b->m ; i++)
    {
        qi = F->Q_perm[i];
        for (j = 0 ; j < b->n ; j++)
        {
            SPEX_MPQ_SET_Z(SPEX_2D(x,  qi, j, mpq),
                                      SPEX_2D(b_new,  i, j, mpz));
            SPEX_MPQ_DIV(SPEX_2D(x,  qi, j, mpq),
                                    SPEX_2D(x,  qi, j, mpq), b_new->scale);
        }
    }

    
    //--------------------------------------------------------------------------
    // Return result and free workspace
    //--------------------------------------------------------------------------
    (*x_handle)=x;//b_new;

    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
