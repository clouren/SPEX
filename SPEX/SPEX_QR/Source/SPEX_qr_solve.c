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

# include "spex_qr_internal.h"
# include "spex_lu_internal.h"

SPEX_info SPEX_qr_solve
(
    SPEX_matrix *x_handle, // Solution
    SPEX_matrix R,        // Upper triangular matrix
    SPEX_matrix Q,
    SPEX_matrix b,        // Q^T * b
    SPEX_matrix rhos,
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

    SPEX_matrix b_new = NULL;
    int64_t k, p, i,j;
    mpz_t det;

    // b->new has Q->n rows and b->n columns
    SPEX_CHECK(SPEX_matrix_allocate(&b_new, SPEX_DENSE, SPEX_MPZ, Q->n, b->n, Q->n*b->n,
        false, true, NULL));
    //det = &(R->x.mpz[R->n-1]); //det = &(F->rhos->x.mpz[F->L->n-1]);
    SPEX_MPZ_INIT(det);
    SPEX_MPZ_SET(det,rhos->x.mpz[R->n-1]);
    // Need to compute b_new[i] = R(n,n)* Q'[i,:] dot b[i]
    // This is equivalent to b_new[i] = R(n,n)* Q[:,i] dot b[i]
    
    // Iterate across every RHS vector
    for (k = 0; k < b->n; k++)
    {
        // Compute b[j,k]
        for(j=0;j<Q->n;j++)
        {
            for(p=Q->p[j]; p < Q->p[j+1]; p++)
            {
                i=Q->i[p];
                SPEX_MPZ_ADDMUL(SPEX_2D(b_new, j, k, mpz),Q->x.mpz[p],SPEX_2D(b, i, k, mpz));
            }
            SPEX_MPZ_MUL (SPEX_2D(b_new, j, k, mpz),SPEX_2D(b_new, j, k, mpz),det);
        }
    }
    
    /*printf("orint b_new sparse:\n");
    SPEX_matrix_check(b_new, option);*/
    //backwards substitution
    //Solves Rx=b_new (overwrites b_new into x)
    SPEX_CHECK (spex_left_lu_back_sub(R,b_new));
    
    (*x_handle)=b_new;
    return SPEX_OK;
}
