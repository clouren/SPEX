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
    SPEX_matrix* x_handle // Solution
    SPEX_matrix R,        // Upper triangular matrix
    SPEX_matrix Q,
    SPEX_matrix b,        // Q^T * b
)
{
    SPEX_info info;
    // Check inputs, the number of columns of Q' must equal the number of rows of b
    // Also, both matrices must be mpz and dense
    ASSERT( Q->m == b->m);
    ASSERT( Q->type == SPEX_MPZ);
    ASSERT( b->type == SPEX_MPZ);
    ASSERT( Q->kind == SPEX_CSC);
    ASSERT( b->kind == SPEX_CSC);

    SPEX_matrix b_new = NULL;
    int64_t k, p, i;
    mpz_t det;

    // b->new has Q->n rows and b->n columns
    SPEX_CHECK(SPEX_matrix_allocate(&b_new, SPEX_CSC, SPEX_MPZ, Q->n, b->n, Q->n*b->n,
        false, true, NULL));
    SPEX_MPZ_SET(det, R->x.mpz[R->nz]);
    
    // Need to compute b_new[i] = R(n,n)* Q'[i,:] dot b[i]
    // This is equivalent to b_new[i] = R(n,n)* Q[:,i] dot b[i]
    
    // Iterate across every RHS vector
    for (k = 0; k < b->n; k++)
    {
        // Compute b[i,k]
        for(p =b->p[k]; p< b->p[k+1];p++)
        {
            i=b->i[p];
            SPEX_CHECK (spex_dot_product(b_new->x.mpz[i],Q, i, b, k,option));
            SPEX_MPZ_PROD (b_new->x.mpz[i],b_new->x.mpz[i],det);
        }
    }
    
    //backwards substitution
    //Solves Rx=b_new (overwrites b_new into x)
    SPEX_CHECK (spex_left_lu_back_sub(R,b_new));
    
    (*x_handle)=b_new;
}
