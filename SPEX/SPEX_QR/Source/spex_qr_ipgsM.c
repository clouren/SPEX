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

int binarySearch(int64_t *arr, int low, int high, int key)
{
    if (high < low)
        return -1;
    int mid = (low + high) / 2; /*low + (high - low)/2;*/
    if (key == arr[mid])
        return mid;
    if (key > arr[mid])
        return binarySearch(arr, (mid + 1), high, key);
    return binarySearch(arr, low, (mid - 1), key);
}


SPEX_info spex_qr_ipgsM
(
    SPEX_matrix R,    
    SPEX_matrix Q,    
    SPEX_matrix rhos,         // sequence of pivots
    int64_t *Qj,
    int64_t *col,
    const int64_t j,          // Row of R to compute (col j+1 of Q will also be computed)
    const SPEX_matrix A,      // Matrix to be factored
    int64_t *h,
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
    
    /*printf("qj for iter %ld\n",j);
    for(i=0;i<m;i++)
    {
        printf("%ld\n",Qj[i]);
    }*/
    
    // Update columns j+1 to n of Q (column j+1 is finalized after this)
    for (k=j+2;k<n;k++)
    {
        //find R(j+1,k)
        pR=binarySearch(R->i, R->p[j], R->p[j+1]-1,k);
        printf("k %ld,pR %ld\n",k,pR);
        if(pR==-1) continue;
        
        //for all elements in column k of Q
        for (pQ =Q->p[k]; pQ < Q->p[k+1]; pQ++)
        {
            iQ=Q->i[pQ];
            prev=Qj[iQ];
            printf("pQ %ld iQ %ld prev %ld\n",pQ,iQ,prev);
            
            //fix logic
            if(prev==-1 && k==j+1)
            {
                printf("last hist %ld %ld %ld %ld\n",pQ, h[pQ],k,j);
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(k>1 && h[pQ]>0)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]]);
                }
                final[iQ]=pQ;
            }
            
            
            
            //check if column j of Q had a nonzero element in row iQ
            if(prev==-1) continue;
            SPEX_MPZ_SGN(&sgn, Q->x.mpz[prev]);
            if(sgn==0) continue;
            
            //history update
            if(k>h[pQ]+1) //"an update of Q(i,j)" has been skipped because R(i,l) is zero or Q(i,l) is zero
            {
                printf("hist %ld %ld\n",pQ, h[pQ]);
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(i>1 && h[pQ%m]>-1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]]);
                }
                
            }
            
            
            //IPGE update
            SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]); //Q(i,k)==0 this can happen when A(i,k)=0, but Q(i,k) is nonzero
            if(sgn==0)
            {
                //Q(i,j)=-R(k,j)*Q(i,k)/rho^(k-1)
                // Q[pQ] = Q[pQ] - R[pR]*Q[prev]
                SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[prev]); 
                if(j>=1) //TO CHECK, is this needed?
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                }
            }
            else
            {
                //printf("j %ld\n",j);
                //Q(i,j)=(rho^k*Q(i,j)-R(k,j)*Q(i,k))/rho^(k-1)
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                // Q[pQ] = Q[pQ] - R[pR]*Q[k]
                SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[prev]);
                if(j>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                }
                
            }
            
            h[pQ]=k;
            
            
            if(k==j+1)//could probablu avoid having this here with some nice rearranging
            {
                final[iQ]=pQ;
            }
            
        }
    }
    //for k=j+1
    k=j+1;
        //find R(j+1,k)
        pR=binarySearch(R->i, R->p[j], R->p[j+1]-1,k);
        //printf("k %ld,pR %ld\n",k,pR);


        //for all elements in column k of Q
        for (pQ =Q->p[k]; pQ < Q->p[k+1]; pQ++)
        {
            iQ=Q->i[pQ];
            prev=Qj[iQ];
            printf("k=j+1 pQ %ld iQ %ld prev %ld\n",pQ,iQ,prev);
            
            //fix logic
            if(prev==-1 || pR==-1)
            {
                //check if column j of Q had a nonzero element in row iQ
            //SPEX_MPZ_SGN(&sgn, Q->x.mpz[prev]);
            //if(sgn==0) continue;
                printf("last hist %ld %ld %ld %ld\n",pQ, h[pQ],k,j);
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(k>1 && h[pQ]>0)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]]);
                }
                final[iQ]=pQ;
            }
            else
            {
            
            if(k>h[pQ]+1) //"an update of Q(i,j)" has been skipped because R(i,l) is zero or Q(i,l) is zero
            {
                printf("hist %ld %ld\n",pQ, h[pQ]);
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(i>1 && h[pQ%m]>-1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]]);
                }
                
            }
            
            if(col[iQ]!=j)
            {
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(j>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                }
                continue;
            }
            
            //IPGE update
            SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]); //Q(i,k)==0 this can happen when A(i,k)=0, but Q(i,k) is nonzero
            if(sgn==0)
            {
                //Q(i,j)=-R(k,j)*Q(i,k)/rho^(k-1)
                // Q[pQ] = Q[pQ] - R[pR]*Q[prev]
                SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[prev]); 
                if(j>=1) //TO CHECK, is this needed?
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                }
            }
            else
            {
                //printf("j %ld\n",j);
                //Q(i,j)=(rho^k*Q(i,j)-R(k,j)*Q(i,k))/rho^(k-1)
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                // Q[pQ] = Q[pQ] - R[pR]*Q[k]
                SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[prev]);
                if(j>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j-1]);
                }
                
            }
            }
            h[pQ]=k;
            
            
            
            Qj[iQ]=pQ;
            col[iQ]=k;
            
            
        }
    
    
   /* for(i=0;i<m;i++)
    {
        Qj[i]=final[i];
    }*/
    //SPEX_matrix_check(Q, option);
    /*for(i=0;i<m;i++)
    {
        printf("%ld \n",vec[i]);
    }
    printf("done ipgs\n");*/
    /*printf("j+1 : %ld\n",j+1); 
    for(i=0;i<Q->nz;i++)
    {
        printf("%ld %ld \n",i,Prev[i]);
    }*/

    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
