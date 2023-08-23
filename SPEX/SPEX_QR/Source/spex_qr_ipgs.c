//------------------------------------------------------------------------------
// SPEX_QR/Source/spex_qr_ipgs.c: Integer Preserving Gram-Schmidt
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* Given the non-zero pattern of the j-th row of R, this code performs one iteration of 
 * REF QR via Integer-preserving Gram-Schmidt to obtain row j of R and column j+1 of Q
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
    SPEX_matrix R,    
    SPEX_matrix Q,    
    SPEX_matrix rhos,         // sequence of pivots
    int64_t *Qj,
    int64_t *col,
    const int64_t j,          // Row of R to compute (col j+1 of Q will also be computed)
    const SPEX_matrix A,      // Matrix to be factored
    int64_t *h,
    bool isZeros,
    SPEX_options option
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
    int64_t p, pQ, pR, iR, top, x,l, prev, iQ,k,i;
    int sgn;
    int64_t *final;
    
    isZeros=true; //start by assuming column of Q is linearly dependent
   
    final = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    
    for(k=0;k<m;k++)
    {
        final[k]=-1;
    }
    // Compute row k of R
    for (pR =R->p[j];pR <R->p[j+1];pR++)
    {
        // Obtain the index of the current nonzero
        i = R->i[pR];//column number where j is row number
        // R(j,i) = Q(:,j) dot A(:,i)
        SPEX_CHECK(spex_dot_product(R->x.mpz[pR],Q, j, A, i, option)); 
    }
    SPEX_MPZ_SET(rhos->x.mpz[j],R->x.mpz[R->p[j]]); //rhos stores the diagonal of R
    //SPEX_matrix_check(rhos,option);

    // Update columns j+1 to n of Q (column j+1 is finalized after this)
    for (pR =R->p[j]; pR < R->p[j+1]; pR++)//the +1 is so that we don't use the diagonal
    {
        i = R->i[pR];
        if(i<=j+1) continue;//when i=j+1 (Qj+i is worked on at the end)
    
        //for all elements in column k of Q
        for (pQ =Q->p[i]; pQ < Q->p[i+1]; pQ++)
        {
            iQ=Q->i[pQ];
            prev=Qj[iQ];

            //check if column j of Q had a nonzero element in row iQ
            if(prev==-1) continue;
            SPEX_MPZ_SGN(&sgn, Q->x.mpz[prev]);
            if(sgn==0) continue;
            if(col[iQ]!=j) continue; 
            
            //history update
            if(j+1>h[pQ]+1) //"an update of Q(i,j)" has been skipped because R(i,l) is zero or Q(i,l) is zero
            {
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                if(h[pQ]>0)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]-1]);
                }
                
            }
 
            //IPGE update
            SPEX_CHECK(spex_ipge(Q,R,rhos,pQ,pR,j-1,j,prev,option));
            
            h[pQ]=j+1;
            
        }
    }
    
    
    //for k=j+1
    //IPGE and finalize the first vector
    k=j+1;
    //find R(j+1,k)
    //pR=binarySearch(R->i, R->p[j], R->p[j+1]-1,k);
    //printf("k %ld,pR %ld\n",k,pR);
    for(pR =R->p[j]; pR < R->p[j+1]; pR++)
    {
        i=R->i[pR];
        //printf("i %ld pR %ld \n",i, pR);
        if(i>=k) break; //should happen in the first couple of elements of R[j]
    }
    //printf("k %ld,pR %ld, i %ld\n",k,pR,i );
    //for all elements in column k of Q
    for (pQ =Q->p[k]; pQ < Q->p[k+1]; pQ++)
    {
        iQ=Q->i[pQ];
        prev=Qj[iQ];
        
        //fix logic
        if(prev==-1 || i!=k)
        {
            //check if column j of Q had a nonzero element in row iQ
            //SPEX_MPZ_SGN(&sgn, Q->x.mpz[prev]);
            //if(sgn==0) continue;
            //Q(i,j)=rho^()*Q(i,k)/rho^()
            // Q[pQ] = x[pQ] * rho[i]
            
            SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
            if(k>1 && h[pQ]>0)
            {
                SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]-1]);
            }
            final[iQ]=pQ;
            h[pQ]=k;
        }
        else
        {
            
            if(j+1>h[pQ]+1) //"an update of Q(i,j)" has been skipped because R(i,l) is zero or Q(i,l) is zero
            {
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                if(h[pQ]>0)
                {
                    
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]-1]);
                }
                
            }
            
            if(col[iQ]!=j)
            {
                printf("HERE\n");
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(j>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                }
                h[pQ]=k;
                continue;
            }
            
            //IPGE update
            SPEX_CHECK(spex_ipge(Q,R,rhos,pQ,pR,j-1,j,prev,option));
        }
        h[pQ]=k;
        
        final[iQ]=pQ;
        col[iQ]=k;
        
        SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]);
        if(sgn!=0)
        {
            isZeros=false;
        }
    }
    
    
    for(i=0;i<m;i++)
    {
        Qj[i]=final[i];
    }
    //SPEX_matrix_check(Q,option);
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
