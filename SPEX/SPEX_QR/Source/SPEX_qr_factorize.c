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
    SPEX_matrix_free(&(RT),option); \
    SPEX_free(h);                   \
    SPEX_free(col);                 \
    SPEX_free(Qk);                   \
    SPEX_free(ldCols);                   \
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
    int64_t n=A->n, m=A->m, k,i,pQ,p, iQ, pR, rankDeficient=0;
    SPEX_factorization F = NULL ;
    SPEX_matrix RT;
    int64_t *h, *Qk, *col;
    bool isZeros=false, *ldCols;
    
    
     // Allocate memory for the factorization
    F = (SPEX_factorization) SPEX_calloc(1, sizeof(SPEX_factorization_struct));
    if (F == NULL) return SPEX_OUT_OF_MEMORY;

    // set factorization kind
    F->kind = SPEX_QR_FACTORIZATION;

    // Allocate and set scale_for_A
    SPEX_MPQ_INIT(F->scale_for_A);
    SPEX_MPQ_SET(F->scale_for_A, A->scale);

    //--------------------------------------------------------------------------
    // Numerically permute matrix A, that is apply the column ordering from
    // the symbolic analysis step to get the permuted matrix AQ.
    //--------------------------------------------------------------------------
    
    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));

    SPEX_CHECK(spex_qr_nonzero_structure(&RT, &F->Q, A, S, option));
    //SPEX_matrix_check(F->Q, option);

    h = (int64_t*) SPEX_calloc((F->Q->nz),sizeof(int64_t));
    Qk = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    col = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    // initialize workspace history array
    /*for (i = 0; i < F->Q->nz; i++)
    {
        h[i] = 0;
    }*/
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
    
    ldCols = (bool*) SPEX_calloc((n),sizeof(bool));//bool?
   
    //--------------------------------------------------------------------------
    // Perform IPGS to get Q and R
    //--------------------------------------------------------------------------
    for (k=0;k<n-1;k++)
    {
        //when the kth column of Q is all zeros (it is linearly dependent)
        //then the kth row of R is all zeros too and you skip operations on k
        if(isZeros)
        {
            printf("here\n");
            ldCols[k]=true;//kth pivot of R is zeros, kth column of Q is ld
            // row k of R is zeros
            for (pR =RT->p[k];pR <RT->p[k+1];pR++)
            {
                SPEX_MPZ_SET_UI(RT->x.mpz[pR],0); 
            }
            // Set the kth pivot to be equal to the k-1th pivot for computations
            SPEX_MPZ_SET(F->rhos->x.mpz[k],F->rhos->x.mpz[k-1]);
            for(i=0;i<m;i++)
            {
                Qk[i]=-1;
            }
            // Finalize Q k+1 (it keeps its previous values)
            for (pQ =F->Q->p[k+1]; pQ < F->Q->p[k+2]; pQ++)
            {
                if(h[pQ]<k)
                {
                    SPEX_CHECK(spex_history_update(F->Q,F->rhos,pQ,k-1,h[pQ],h[pQ]-1,0,option)); //TODO check
                }
                else
                {
                    h[pQ]=k+1;
                }
                iQ=F->Q->i[pQ];
                Qk[iQ]=pQ;
                col[iQ]=k+1;
            }
            // Update the entire history vector
            for(i = pQ; i < F->Q->nz; i++)
            {
                if(h[i]==k)
                {
                    h[i]=k+1;
                }
            }
            rankDeficient++;
            
            isZeros=false;
        }
        else
        {
            // Integer-preserving Gram-Schmidt
            SPEX_CHECK(spex_qr_ipgs(RT, F->Q, F->rhos, Qk,col, k, A, h,
                                     &isZeros, S->Q_perm, option));
        }
        //SPEX_matrix_check(F->Q, option);
        //SPEX_matrix_check(RT, option);
    }
    
    // Finish R (get the last element/pivot)
    if(isZeros)
    {
        ldCols[k]=true;
        SPEX_MPZ_SET_UI(RT->x.mpz[RT->p[n]-1],0); 
        SPEX_MPZ_SET(F->rhos->x.mpz[n-1],F->rhos->x.mpz[n-2]);
    }
    else
    {
        SPEX_CHECK(spex_dot_product(RT->x.mpz[RT->p[n]-1],F->Q, n-1, A, S->Q_perm[n-1], option)); 
        SPEX_MPZ_SET(F->rhos->x.mpz[n-1],RT->x.mpz[RT->p[n]-1]);    
    }
    //--------------------------------------------------------------------------
    // Get rank revealing permutation 
    //--------------------------------------------------------------------------
    if(rankDeficient>0)
    {
        int64_t *Pi_perm;
        int64_t iZero=n-1, iNon=0;
        
        SPEX_matrix QPi=NULL, RTPi=NULL;
        
        F->rank=n-rankDeficient;
        printf("rank %ld %ld\n",F->rank,rankDeficient);
        
        Pi_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) ); //if we don't permute A we don't need Pi perm unless we want to store both F->pi and F->Q
        F->Q_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        if (!(F->Q_perm))
        {
            // out of memory: free everything and return
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }

        for(k=0;k<n;k++)
        {
            if(ldCols[k]) //ldCols[k] is true when the k col is linearly dependent
            {
                //Pi_perm[iZero]=k;
                F->Q_perm[iZero]=S->Q_perm[k];
                iZero--;
            }
            else
            {
                //Pi_perm[iNon]=k;
                F->Q_perm[iNon]=S->Q_perm[k];
                iNon++;
            }
        }
        
        /*SPEX_CHECK( spex_qr_permute_A(&QPi, F->Q, true, Pi_perm, option) );
        SPEX_CHECK( spex_qr_permute_A(&RTPi, RT, true, Pi_perm, option) );
        
        SPEX_CHECK(SPEX_transpose(&F->R,RTPi,option));
        F->R->nz=RT->p[n]-1;
        
        F->Q=QPi;//TODO check memleaks, maybe should coypy
        */
              
    }
    else
    {
        
        // column permutation, to be copied from S->Q_perm
        F->Q_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        if (!(F->Q_perm))
        {
            // out of memory: free everything and return
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }

        // Copy column permutation from symbolic analysis to factorization
        memcpy(F->Q_perm, S->Q_perm, n*sizeof(int64_t));

    }
    
    
    SPEX_CHECK(SPEX_transpose(&F->R,RT,option));
    F->R->nz=RT->p[n]-1;

    //--------------------------------------------------------------------------
    // Return result and free workspace
    //--------------------------------------------------------------------------

    (*F_handle)=F;
    
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
