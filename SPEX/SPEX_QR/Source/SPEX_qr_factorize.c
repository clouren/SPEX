//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_qr_factorize.c: QR factorization
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* Purpose: This function performs the REF QR factorization via 
 * Integer-preserving Gram-Schmidt
 *
 * Input arguments of the function:
 *
 * F_handle:    Handle to the factorization struct. Null on input.
 *              On output, contains a pointer to the factorization.
 *
 * A:           User's input matrix. Must be SPEX_MPZ and SPEX_CSC.
 *
 * S:           Symbolic analysis struct for QR factorization.
 *              On input it contains the column elimination tree and
 *              the number of nonzeros in R.
 *
 * option:      Command options.
 */

# include "spex_qr_internal.h"

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_matrix_free(&(RT),option); \
    SPEX_free(h);                   \
    SPEX_free(Qk);                  \
    SPEX_free(ldCols);              \
}

# define SPEX_FREE_ALL                \
{                                     \
    SPEX_FREE_WORKSPACE               \
    SPEX_factorization_free(&F, NULL);\
}



SPEX_info SPEX_qr_factorize
(
    // Output
    SPEX_factorization *F_handle,   // QR factorization struct
    //Input
    const SPEX_matrix A,            // Matrix to be factored. Must be SPEX_MPZ
                                    // and SPEX_CSC
    SPEX_symbolic_analysis S,       // Symbolic analysis struct containing the
                                    // column elimination tree of A, the column
                                    // permutation, and number of nonzeros in R
    const SPEX_options option       // command options.
)
{
    SPEX_info info;

    // Declare variables
    int64_t n=A->n, m=A->m, k, i, pQ, p, iQ, pR;
    SPEX_factorization F = NULL ;
    SPEX_matrix RT, Q;
    int64_t *h, *Qk;
    // Varibles needed to compute the rank of a matrix.
    // assume matrix is full rank. isZeros is true if a column is linearly dependent
    // ldCols keeps track of linearly dependent columns
    bool isZeros=false, *ldCols;
    int64_t rank=n;
    
    
     // Allocate memory for the factorization
    F = (SPEX_factorization) SPEX_calloc(1, sizeof(SPEX_factorization_struct));
    if (F == NULL) return SPEX_OUT_OF_MEMORY;

    // set factorization kind
    F->kind = SPEX_QR_FACTORIZATION;

    // Allocate and set scale_for_A
    SPEX_MPQ_INIT(F->scale_for_A);
    SPEX_MPQ_SET(F->scale_for_A, A->scale);

    //--------------------------------------------------------------------------
    // Allocate and compute the nonzero structure of Q and R
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_qr_nonzero_structure(&RT, &Q, A, S, option));
    
    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));

    //--------------------------------------------------------------------------
    // Allocate and initialize supporting vectors
    //--------------------------------------------------------------------------
    //TODO add explanation for Qk
    h = (int64_t*) SPEX_calloc((Q->nz),sizeof(int64_t)); //history matrix
    Qk = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    for(k=0;k<m;k++)
    {
        Qk[k]=-1;
    }
    for(k=Q->p[0];k<Q->p[1];k++)
    {
        Qk[Q->i[k]]=k;
    }
    
    ldCols = (bool*) SPEX_calloc((n),sizeof(bool));
   
    //--------------------------------------------------------------------------
    // Perform IPGS to get Q and R
    //--------------------------------------------------------------------------
    for (k=0;k<n-1;k++)
    {
        //when the kth column of Q is all zeros (it is linearly dependent)
        //then the kth row of R is all zeros too and you skip operations on k
        if(isZeros)
        {
            ldCols[k]=true;//kth pivot of R is zeros, kth column of Q is ld
            // row k of R is zeros
            for (pR =RT->p[k];pR <RT->p[k+1];pR++)
            {
                SPEX_MPZ_SET_UI(RT->x.mpz[pR],0); 
            }
            // Set the kth pivot to be equal to the k-1th pivot for computations
            SPEX_MPZ_SET(F->rhos->x.mpz[k],F->rhos->x.mpz[k-1]);
            for(i=0;i<m;i++) //TODO change to For (I = A->p[k-2]; I < A->p[k-1]; I++) "erase prev nonzeros"
            {
                Qk[i]=-1;
            }
            // Finalize Q k+1 (it keeps its previous values)
            for (pQ =Q->p[k+1]; pQ < Q->p[k+2]; pQ++)
            {
                //History update
                if(h[pQ]<k)
                {
                    SPEX_CHECK(spex_history_update(Q,F->rhos,pQ,k-1,h[pQ],
                                                    h[pQ]-1,0,option));
                }
                else
                {
                    h[pQ]=k+1;
                }
                iQ=Q->i[pQ];
                Qk[iQ]=pQ;
            }
            // Update the history vector
            for(i = pQ; i < Q->nz; i++)
            {
                if(h[i]==k)
                {
                    h[i]=k+1;
                }
            }

            rank--;
            
            isZeros=false;
        }
        else
        {
            // Integer-preserving Gram-Schmidt
            SPEX_CHECK(spex_qr_ipgs(RT, Q, F->rhos, Qk, h, &isZeros, k, A,
                                     S->Q_perm, option));
        }

    }
    
    // Finalize R (get the last element/pivot)
    if(isZeros)
    {
        ldCols[k]=true;
        SPEX_MPZ_SET_UI(RT->x.mpz[RT->p[n]-1],0); 
        SPEX_MPZ_SET(F->rhos->x.mpz[n-1],F->rhos->x.mpz[n-2]);
    }
    else
    {
        SPEX_CHECK(spex_dot_product(RT->x.mpz[RT->p[n]-1],Q, n-1, A, 
                                       S->Q_perm[n-1], option)); 
        SPEX_MPZ_SET(F->rhos->x.mpz[n-1],RT->x.mpz[RT->p[n]-1]);    
    }

    //--------------------------------------------------------------------------
    // Get rank revealing permutation 
    //--------------------------------------------------------------------------
    if(rank!=n)
    {
        int64_t *Pi_perm;
        int64_t iZero=n-1, iNon=0;
        
        SPEX_matrix RTPi=NULL;
        
        F->rank=rank;
        
        Pi_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        F->Q_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        if (!(F->Q_perm))
        {
            // out of memory: free everything and return
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }

        // TODO coment more, explain iZero, iNon, change to iLD iLI
        for(k=0;k<n;k++)
        {
            if(ldCols[k]) //ldCols[k] is true when the k col is linearly dependent
            {
                Pi_perm[iZero]=k;
                F->Q_perm[iZero]=S->Q_perm[k];
                iZero--;
            }
            else
            {
                Pi_perm[iNon]=k;
                F->Q_perm[iNon]=S->Q_perm[k];
                iNon++;
            }
        }

        // Permute Q and RT
        // Zero columns of Q will be on the right 
        // Zero rows of R will be on the bottom (Zero columns of RT will be on the right)
        SPEX_CHECK( spex_qr_permute_A(&F->Q, Q, true, Pi_perm, option) );
        SPEX_CHECK( spex_qr_permute_A(&RTPi, RT, true, Pi_perm, option) );
        
        
        SPEX_CHECK(SPEX_transpose(&F->R,RTPi,option));
        F->R->nz=RT->p[n]-1;
        
        SPEX_matrix_free(&Q,option);
    }
    else
    {
        F->rank=n; //matrix has full rank
        
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
        
        F->Q=Q;
        SPEX_CHECK(SPEX_transpose(&F->R,RT,option));
        F->R->nz=RT->p[n]-1;

    }
    

    //--------------------------------------------------------------------------
    // Return result and free workspace
    //--------------------------------------------------------------------------

    (*F_handle)=F;
    
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
