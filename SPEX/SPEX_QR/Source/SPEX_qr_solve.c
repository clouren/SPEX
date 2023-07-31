//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_qr_solve.c: Solve exactly x = R \ Q^T b
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs the sparse R \ Q^T b. First it computes R(n,n)*Q^T b has already
 * been computed. Then does
 * back solve.  Returns x as integer.
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

    SPEX_matrix b_new = NULL, b_perm=NULL;
    int64_t k, p, i,j;
    mpz_t det;

    // b->new has Q->n rows and b->n columns
    SPEX_CHECK(SPEX_matrix_allocate(&b_new, SPEX_DENSE, SPEX_MPZ, F->Q->n, b->n, F->Q->n*b->n,
        false, true, NULL));
    SPEX_MPZ_INIT(det);
    SPEX_MPZ_SET(det,F->rhos->x.mpz[F->R->n-1]);
    // Need to compute b_new[i] = R(n,n)* Q'[i,:] dot b[i]
    // This is equivalent to b_new[i] = R(n,n)* Q[:,i] dot b[i]
    
    SPEX_CHECK (spex_permute_dense_matrix (&b_perm, b, F->Q_perm, option)); 
    // Iterate across every RHS vector
    for (k = 0; k < b->n; k++)
    {
        // Compute b[j,k]
        for(j=0;j<F->Q->n;j++)
        {
            for(p=F->Q->p[j]; p < F->Q->p[j+1]; p++) //dot product probably not good sparse FIXME
            {
                i=F->Q->i[p];
                SPEX_MPZ_ADDMUL(SPEX_2D(b_new, j, k, mpz),F->Q->x.mpz[p],SPEX_2D(b_perm, i, k, mpz));
            }
            SPEX_MPZ_MUL (SPEX_2D(b_new, j, k, mpz),SPEX_2D(b_new, j, k, mpz),det);
        }
    }
    
    //backwards substitution
    //Solves Rx=b_new (overwrites b_new into x)
    SPEX_CHECK (spex_left_lu_back_sub(F->R,b_new));
    
    (*x_handle)=b_new;
    return SPEX_OK;
}
