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

    if (!spex_initialized())
    {
        return SPEX_PANIC;
    }

    // Check inputs for NULL
    if (!F_handle || !A || !S)
    {
        return (SPEX_INCORRECT_INPUT);
    }

    // Ensure inputs are in the correct format
    if (A->kind != SPEX_CSC || A->type != SPEX_MPZ
        || S->kind != SPEX_QR_FACTORIZATION)
    {
        return (SPEX_INCORRECT_INPUT);
    }

    SPEX_factorization_algorithm algo = SPEX_OPTION_ALGORITHM(option);
    if(algo!=SPEX_QR_GS)
    {
        return SPEX_INCORRECT_ALGORITHM;
    }

    // Declare variables
    int64_t n=A->n, m=A->m, k, i, pQ, p, iQ, pR;
    SPEX_factorization F = NULL ;
    SPEX_matrix RT, Q;
    int64_t *h, *Qk;
    // Varibles needed to compute the rank of a matrix.
    // assume matrix is full rank. isZeros is true if a column is linearly 
    // dependent, ldCols keeps track of linearly dependent columns
    bool isZeros=false, *ldCols;
    int64_t rank=n;
    int sgn;
    
    clock_t start, end;
    double times;
    
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
 
    /*for (i = 0; i < n+1; ++i)
    {
        printf("i %ld Q %ld R %ld\n",i,Q->p[i],RT->p[i]);
    }*/
    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));

    //--------------------------------------------------------------------------
    // Allocate and initialize supporting vectors
    //--------------------------------------------------------------------------
    h = (int64_t*) SPEX_calloc((Q->nz),sizeof(int64_t)); //history matrix
    // Qk contains the position of the nonzero elements in each row of the 
    // latest column of Q to be finalized (or -1 if the corresponding element 
    // in that row is symbolically zero
    Qk = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    for(k=0;k<m;k++)
    {
        Qk[k]=-1;
    }
    for(k=Q->p[0];k<Q->p[1];k++) // Q(:,0)=A(:,0) first column of Q is finalized
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
            //printf("%ld in isZeros\n",k);
            ldCols[k]=true;//kth pivot of R is zeros, kth column of Q is ld

            isZeros=true;

            // row k of R is zeros
            for (pR =RT->p[k];pR <RT->p[k+1];pR++)
            {
                SPEX_MPZ_SET_UI(RT->x.mpz[pR],0); 
            }
            // Set the kth pivot to be equal to the k-1th pivot for computations
            SPEX_MPZ_SET(F->rhos->x.mpz[k],F->rhos->x.mpz[k-1]);
            for(i=Q->p[k];i<Q->p[k+1];i++) //erase prev nonzeros
            {
                Qk[Q->i[i]]=-1;
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

                SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]);
                if(sgn!=0)
                {
                    //printf("k %ld\n",k);
                    isZeros=false;
                }
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

            if(k==11){
                //option->print_level = 3;
                //SPEX_matrix_check(F->rhos, option);
            }
        }
        else
        {
            // Integer-preserving Gram-Schmidt
            SPEX_CHECK(spex_qr_ipgs(RT, Q, F->rhos, Qk, h, &isZeros, k, A,
                                     S->Q_perm, option));
            //option->print_level = 3;
            //SPEX_matrix_check(Q, option);
            //SPEX_matrix_check(RT, option);
        }

    }

    
    // Finalize R (get the last element/pivot)
    if(isZeros)
    {
        //printf("last isZeros\n");
        ldCols[k]=true;
        rank--;
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
        // If A is rank deficient then Q and RT have at least one 0 column
        // Here those columns will be permuted to the right side of the 
        // respective matrix. And as such the column permutation of A (Q_perm)
        // will be updated.
        
        int64_t *Pi_perm; // Column permutation for rank deficient matrices
        int64_t *Piinv_perm; // Inverse row permutation for rank deficient matrices
        // Indices for modifying Q_perm (and creating Pi_perm) according to 
        // whether a column of Q is linearly dependent or linearly independent 
        // of the previous columns
        int64_t iLD=n-1, iLI=0, index;
        
        SPEX_matrix RTPi=NULL;
        SPEX_matrix RPi=NULL;
        
        F->rank=rank;
        
        Pi_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        Piinv_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        F->Q_perm = (int64_t*) SPEX_malloc ( n*sizeof(int64_t) );
        if (!(F->Q_perm))
        {
            // out of memory: free everything and return
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }

        // If the k-th column of Q is linearly dependent (ldCols[k]=true), then
        // k needs to be towards the end of Pi_perm and the original Q_perm[k]
        // needs to be towards the end of the updated column permutation
        for(k=0;k<n;k++)
        {
            if(ldCols[k]) //ldCols[k] is true when the k col is linearly dependent
            {
                Pi_perm[iLD]=k;
                //Piinv_perm[n+iLD]=k;
                F->Q_perm[iLD]=S->Q_perm[k];
                iLD--;
            }
            else
            {
                Pi_perm[iLI]=k;
                //Piinv_perm[n-iLI]=k;
                F->Q_perm[iLI]=S->Q_perm[k];
                iLI++;
            }
        }
        
        // Populate pinv
        for (k = 0; k < n; k++)
        {
            index = Pi_perm[k];
            Piinv_perm[index] = k;
        }

        // Permute Q and RT
        // Zero columns of Q will be on the right 
        // Zero rows of R will be on the bottom (Zero columns of RT will be on 
        // the right)
        SPEX_CHECK( spex_qr_permute_A(&F->Q, Q, true, Pi_perm, option) );
        SPEX_CHECK( spex_qr_permute_A2(&RTPi, RT, true, Pi_perm, Piinv_perm, option) );
        SPEX_CHECK(SPEX_transpose(&F->R,RTPi, true, option));
        F->R->nz=RT->p[n]-1;
        
        SPEX_matrix_free(&Q,option);
        SPEX_matrix_free(&RTPi,option);
        SPEX_free(Pi_perm);
        SPEX_free(Piinv_perm);
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
        SPEX_CHECK(SPEX_transpose(&F->R,RT, true, option));
        F->R->nz=RT->p[n]-1;

    }
    

    //--------------------------------------------------------------------------
    // Return result and free workspace
    //--------------------------------------------------------------------------
    (*F_handle)=F;
    
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
