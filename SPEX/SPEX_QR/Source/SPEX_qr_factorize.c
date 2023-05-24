//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_qr_factorize.c: QR factorization
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs the REF QR factorization via Integer-preserving Gram-Schmidt
 */

# include "spex_qr_internal.h"


SPEX_info SPEX_qr_factorize
(
    SPEX_matrix *R_handle,    // Null on input, contains R on output
    SPEX_matrix *Q_handle,    // Null on input, contains Q on output
    SPEX_matrix *rhos_handle,
    const SPEX_matrix A,      // Matrix to be factored
    const SPEX_symbolic_analysis S,
    SPEX_options option
)
{
    SPEX_info info;

    // Declare variables
    int64_t *c, *h, *xi;
    int64_t n=A->n, m=A->m, k,i,pQ;
    SPEX_matrix rhos;
    SPEX_matrix R, Q;

    //Allocate variables

/*    c = (int64_t*) SPEX_malloc(n*sizeof(int64_t));

    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed into the triangular
    // solve algorithm
    h = (int64_t*) SPEX_malloc(n*sizeof(int64_t));

    // xi serves as a global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L
    // for the triangular solve.
    xi = (int64_t*) SPEX_malloc(2*n*sizeof(int64_t));

    if (!h || !xi || !c)
    {
        SPEX_FREE_WORKSPACE;
        return SPEX_OUT_OF_MEMORY;
    }

    // initialize workspace history array
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
    }*/

    SPEX_CHECK (SPEX_matrix_allocate(&(rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));

    /*SPEX_CHECK (SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, n*n, false, false, NULL));
    // Set the column pointers of R
    for (k = 0; k < n; k++)
    {   
        R->p[k] = c[k] = S->cp[k+1];
    }*/
    // we pre-allocate R
    // by performing a symbolic version of the factorization and obtaining the
    // exact nonzero pattern of R

    SPEX_CHECK(spex_qr_pre_factor(&R, A, S));

    /*SPEX_CHECK (SPEX_matrix_allocate(&Q, SPEX_CSC, SPEX_MPZ, m, n, m*n, false, false, NULL));
    //Set col pointers Q THIS IS DENSE
    for (k = 0; k < n; k++)
    {   
        Q->p[k] = k*n;
    }
    // Set Q(:,0)=A(:,0)
    for (pQ =Q->p[0]; pQ < Q->p[1]; pQ++) //if we had a pattern for Q_j this is where it would go
    {
        Q->i[pQ]= pQ % m;
        if(A->i[pQ]==Q->i[pQ])
        {
            SPEX_MPZ_SET(Q->x.mpz[pQ],A->x.mpz[pQ]); 
        }
        else
        {//I don't like this
            SPEX_MPZ_SET_UI(Q->x.mpz[pQ],0);
        }
    }*/
    SPEX_CHECK(SPEX_matrix_copy(&Q, SPEX_CSC, SPEX_MPZ, A, NULL));//this is not the right way to do this because Q will be denser than A
    // Perform IPGS to get Q and R

    for (k=0;k<n-1;k++)
    {
        SPEX_CHECK(spex_qr_ipgs(R, Q, rhos, k, A,/* h, xi, c, S->parent,*/ option));
    }

    //finish R
    SPEX_CHECK(spex_dot_product(R->x.mpz[R->p[n]-1],Q, n-1, A, n-1, option)); 
    SPEX_MPZ_SET(rhos->x.mpz[n-1],R->x.mpz[R->p[n]-1]);

    (*R_handle)=R;
    (*Q_handle)=Q;
    (*rhos_handle)=rhos;
    return SPEX_OK;
}
