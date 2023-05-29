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
    SPEX_factorization *F_handle,
    const SPEX_matrix A,      // Matrix to be factored
    const SPEX_symbolic_analysis S,
    SPEX_options option
)
{
    SPEX_info info;

    // Declare variables
    int64_t *c, *h, *xi;
    int64_t n=A->n, m=A->m, k,i,pQ;
    SPEX_factorization F = NULL ;

    //Allocate variables
    
     // Allocate memory for the factorization
    F = (SPEX_factorization) SPEX_calloc(1, sizeof(SPEX_factorization_struct));
    if (F == NULL) return SPEX_OUT_OF_MEMORY;

    // set factorization kind
    F->kind = SPEX_QR_FACTORIZATION;

    // Allocate and set scale_for_A
    SPEX_MPQ_INIT(F->scale_for_A);
    SPEX_MPQ_SET(F->scale_for_A, A->scale);

    // Inverse pivot ordering
 /*   F->Pinv_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
    // row/column permutation, to be copied from S->P_perm
    F->Q_perm =    (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
    if (!(F->Pinv_perm) || !(F->Q_perm))
    {
        // out of memory: free everything and return
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
/*
    // Copy row/column permutation from symbolic analysis to factorization
    memcpy(F->Q_perm, S->Q_perm, n*sizeof(int64_t));
    memcpy(F->Pinv_perm, S->Pinv_perm, n*sizeof(int64_t));

*/
    //also do estimate here

    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));


    SPEX_CHECK(spex_qr_pre_factor(&F->R, A, S));
    SPEX_matrix_check(F->R, option);
    //SPEX_CHECK(SPEX_matrix_copy(&F->Q, SPEX_CSC, SPEX_MPZ, A, NULL));//this is not the right way to do this because Q will be denser than A
    SPEX_CHECK(spex_qr_pre_Q(&F->Q,A,option));
    // Perform IPGS to get Q and R
   // printf("Q=A:\n");
    //SPEX_matrix_check(F->Q, option);///there is something off with Q that causes a weird free after printing

    for (k=0;k<n-1;k++)
    {
        printf("iteration: %ld\n",k);
        SPEX_CHECK(spex_qr_ipgs(F->R, F->Q, F->rhos, k, A, option));
        
        //here i need to actually assign R(k,:) and Q(:,k+1) with appropiate sizes
    }

    //finish R
    SPEX_CHECK(spex_dot_product(F->R->x.mpz[F->R->p[n]-1],F->Q, n-1, A, n-1, option)); 
    SPEX_MPZ_SET(F->rhos->x.mpz[n-1],F->R->x.mpz[F->R->p[n]-1]);

    (*F_handle)=F;
    return SPEX_OK;
}
