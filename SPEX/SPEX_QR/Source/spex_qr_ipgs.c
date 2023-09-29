//------------------------------------------------------------------------------
// SPEX_QR/Source/spex_qr_ipgs.c: Integer Preserving Gram-Schmidt
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This function obtains row j of R (using a dot product) and column j+1 pf Q 
 * (using Integer-preserving Gram-Schmidt).
 * 
 * 
 * Input/output arguments:
 *
 * R:        Right triangular matrix. On input contains j-1 columns computed
 *           On output contains j columns computed.
 *  
 * Q:        Pair-wise orthogonal matrix. On input contains j columns finalized
 *           and n-j columns in progress. On output contains j+1 columns finalized
 *           and n-j+1 columns in progress.
 *    
 * rhos:     Sequence of pivots. R(k,k)
 *
 * Qj:       Vector containing the pointers to the elements of the j-th column of Q
 *
 * j:        Row of R to compute (col j+1 of Q will be finalized)
 *
 * A:        Matrix to be factored
 *
 * h:        History vector
 * isZeros:  boolean value. On output it's true if column j+1 is linearly dependent
 *           of columns 0 to j of Q, false otherwise. The numerical values of 
 *           linearly dependent columns are zero.
 *
 * Q_perm:   Column permutation
 *
 * option:   Command options
 */

# include "spex_qr_internal.h"

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_FREE(final);                   \
}

# define SPEX_FREE_ALL               \
{                                    \
    SPEX_FREE_WORKSPACE              \
    SPEX_matrix_free(&R, NULL);      \
    SPEX_matrix_free(&Q, NULL);      \
    SPEX_matrix_free(&rhos, NULL);      \
}


SPEX_info spex_qr_ipgs
(
    //Input/Output
    SPEX_matrix R,       // Right triangular matrix
    SPEX_matrix Q,       // Pair-wise orthogonal matrix
    SPEX_matrix rhos,    // sequence of pivots
    int64_t *Qj,         // pointers to elements of the jth column of Q
    int64_t *h,          // History vector
    //Output
    bool *isZeros,       // True if j+1th column of Q is linearly dependent
    //Input
    const int64_t j,     // Row of R to compute (col j+1 of Q will be finalized)
    const SPEX_matrix A, // Matrix to be factored
    const int64_t *Q_perm,     // Column permutation
    const SPEX_options option  // Command options
)
{
    SPEX_info info;
    int64_t m = A->m, n = A->n;
    ASSERT( m >= n); // A should be transposed if not true
    if (m < n)
        return SPEX_PANIC;
    ASSERT( A != NULL);
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_CSC);
    ASSERT (R != NULL);
    ASSERT (Q != NULL);
    
    // Declare variables
    int64_t p, pQ, pR, iR, top, x,l, prev, iQ, k, i;
    int sgn;
    int64_t *final;

    clock_t start, end;
    double times;
    size_t size;
    
    *isZeros=true; //start by assuming column of Q is linearly dependent
   
    final = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));

    //--------------------------------------------------------------------------
    // Compute row j of R, store as column
    //--------------------------------------------------------------------------
    //start = clock();
    for (pR =R->p[j];pR <R->p[j+1];pR++)
    {
        // Obtain the index of the current nonzero
        i = R->i[pR];//column number where j is row number
        // R(j,i) = Q(:,j) dot AQ(:,i)
        SPEX_CHECK(spex_dot_product(R->x.mpz[pR], Q, j, A, Q_perm[i], option)); 
    }
   /* end = clock();
    times=(double) (end - start) / CLOCKS_PER_SEC;
    printf("k  %ld timeR %f",j, times);*/
    //rhos stores the diagonal of R (pivots)
    SPEX_MPZ_SET(rhos->x.mpz[j],R->x.mpz[R->p[j]]);
    
    //--------------------------------------------------------------------------
    // IPGE and finalize column j+1 of Q
    //--------------------------------------------------------------------------
    k=j+1;
    //start=clock();
    // Find the necessary element of R
    //printf("k: %ld\n",k);
    for(pR = R->p[j]; pR < R->p[j+1]; pR++)
    {
        i=R->i[pR];
        if(i>=k) break; //should happen in the first couple of elements of R[j]
    }

    // For all elements in column k of Q
    for (pQ = Q->p[k]; pQ < Q->p[k+1]; pQ++)
    {
        iQ=Q->i[pQ];
        prev=Qj[iQ];
        
        if(prev==-1 || i!=k)
        {
            SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
            if(k>1 && h[pQ]>0)
            {
                SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], 
                                    rhos->x.mpz[h[pQ]-1]);
            }
        }
        else
        {     
            if(j+1>h[pQ]+1)
            {
                //"an update of Q(iQ,j+1)" has been skipped because R(j,i) is zero 
                // or Q(iQ,j) is zero
                SPEX_CHECK(spex_history_update(Q,rhos,pQ,j-1,h[pQ],h[pQ]-1,0,option));
            }
            
            //IPGE update
            SPEX_CHECK(spex_ipge_update(Q,R,rhos,pQ,pR,j-1,j,prev,option));
        }

        // Update history vector
        h[pQ]=k;
        
        // Update the final vectors needed for the next iteration
        final[iQ]=pQ;
        
        // Checks if the entire column k of Q is zeros for rank revealing QR
        SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]);
        if(sgn!=0)
        {
            *isZeros=false;
        }

        size = mpz_sizeinbase(Q->x.mpz[pQ],10);
        //printf("%zu, ",size);
    }
    //printf("\n");
    /*end = clock();
    times=(double) (end - start) / CLOCKS_PER_SEC;
    printf(" timeQj+1 %f",times);

    start=clock();*/
    //--------------------------------------------------------------------------
    // Update columns j+2 to n of Q
    //--------------------------------------------------------------------------
    for (pR =R->p[j]; pR < R->p[j+1]; pR++)
    {
        i = R->i[pR];
        if(i<=j+1) continue;//the j+1 column has already been updated/finalized
    
        //for all elements in column i of Q
        for (pQ =Q->p[i]; pQ < Q->p[i+1]; pQ++)
        {
            iQ=Q->i[pQ];
            prev=Qj[iQ];

            //check if column j of Q had a zero element in row iQ
            if(prev==-1) continue;//simbolic zero
            SPEX_MPZ_SGN(&sgn, Q->x.mpz[prev]);
            if(sgn==0) continue;//numeric zero
        
            //history update
            //"an update of Q(iQ,i)" has been skipped because R(j,i) is zero 
            // or Q(iQ,i-1) is zero
            if(j+1>h[pQ]+1)
            {
                SPEX_CHECK(spex_history_update(Q,rhos,pQ,j-1,h[pQ],h[pQ]-1,0,option));
            }
 
            // IPGE update
            SPEX_CHECK(spex_ipge_update(Q,R,rhos,pQ,pR,j-1,j,prev,option));
            
            // Record changes in history vector
            h[pQ]=j+1;
            /*if(j+1<=5)
            {
                printf("iR %ld iQ %ld pQ %ld h[pQ] %ld\n",i,iQ,pQ,h[pQ]);
            }*/
        }
    }
   /* end = clock();
    times=(double) (end - start) / CLOCKS_PER_SEC;
    printf(" timeQall %f\n",times);*/
    
    // Update the final and col vectors needed for the next iteration
    for(i = Q->p[j]; i < Q->p[j+1]; i++)
    {
        Qj[Q->i[i]]=-1;
    }
    for(i = Q->p[j+1]; i < Q->p[j+2]; i++)
    {
        Qj[Q->i[i]]=final[Q->i[i]];
    }
    //SPEX_matrix_check(Q, option); 
    //--------------------------------------------------------------------------
    // Free workspace
    //--------------------------------------------------------------------------
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
