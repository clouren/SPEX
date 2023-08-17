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

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_matrix_free(&PAQ, option); \
    SPEX_matrix_free(&(F->L),option); \
    SPEX_free(h);                   \
    SPEX_free(col);                 \
    SPEX_free(Qk);                   \
}

# define SPEX_FREE_ALL               \
{                                    \
    SPEX_FREE_WORKSPACE              \
    SPEX_factorization_free(&F, NULL);      \
}



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
    int64_t n=A->n, m=A->m, k,i,pQ,p;
    SPEX_factorization F = NULL ;
    SPEX_matrix PAQ;

    //Allocate variables
    int64_t *h, *Qk, *col;
    
     // Allocate memory for the factorization
    F = (SPEX_factorization) SPEX_calloc(1, sizeof(SPEX_factorization_struct));
    if (F == NULL) return SPEX_OUT_OF_MEMORY;

    // set factorization kind
    F->kind = SPEX_QR_FACTORIZATION;

    // Allocate and set scale_for_A
    SPEX_MPQ_INIT(F->scale_for_A);
    SPEX_MPQ_SET(F->scale_for_A, A->scale);

    // Inverse pivot ordering
    F->Pinv_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) ); 
    // row/column permutation, to be copied from S->P_perm
    F->Q_perm =    (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
    if (!(F->Pinv_perm) || !(F->Q_perm))
    {
        // out of memory: free everything and return
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // Copy row/column permutation from symbolic analysis to factorization
    memcpy(F->Q_perm, S->Q_perm, n*sizeof(int64_t));
    memcpy(F->Pinv_perm, S->Pinv_perm, n*sizeof(int64_t));

    //--------------------------------------------------------------------------
    // Numerically permute matrix A, that is apply the row/column ordering from
    // the symbolic analysis step to get the permuted matrix PAP.
    //--------------------------------------------------------------------------

    SPEX_CHECK( spex_qr_permute_A(&PAQ, A, true, S, option) );

    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));


    SPEX_CHECK(spex_qr_nonzero_structure(&F->L, &F->Q, PAQ, S));
    
    //SPEX_matrix_check(F->Q, option);

    h = (int64_t*) SPEX_malloc((F->Q->nz)*sizeof(int64_t));
    Qk = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    col = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    // initialize workspace history array
    for (i = 0; i < F->Q->nz; i++)
    {
        h[i] = 0;
    }
    for(k=0;k<m;k++)
    {
        Qk[k]=-1;
        col[k]=-1;
    }
    for(k=F->Q->p[0];k<F->Q->p[1];k++)
    {
        Qk[F->Q->i[k]]=k;
        col[F->Q->i[k]]=0;
    }
    
    
    // Perform IPGS to get Q and R
    for (k=0;k<n-1;k++)
    {
        SPEX_CHECK(spex_qr_ipgs(F->L, F->Q, F->rhos, Qk,col, k, PAQ, h, option));
    }

    //finish R
    /*SPEX_MPZ_INIT(F->R->x.mpz[F->R->p[n]-1]);
    SPEX_CHECK(spex_dot_product(F->R->x.mpz[F->R->p[n]-1],F->Q, n-1, PAQ, n-1, option)); 
    SPEX_MPZ_SET(F->rhos->x.mpz[n-1],F->R->x.mpz[F->R->p[n]-1]);
    F->R->nz=F->R->p[n]-1;*/
    SPEX_MPZ_INIT(F->L->x.mpz[F->L->p[n]-1]);
    SPEX_CHECK(spex_dot_product(F->L->x.mpz[F->L->p[n]-1],F->Q, n-1, PAQ, n-1, option)); 
    SPEX_MPZ_SET(F->rhos->x.mpz[n-1],F->L->x.mpz[F->L->p[n]-1]);
    F->L->nz=F->L->p[n]-1;
    
    
    SPEX_CHECK(SPEX_transpose(&F->R,F->L,option));
    F->R->nz=F->L->p[n]-1;


    (*F_handle)=F;
    
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
