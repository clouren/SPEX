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
    SPEX_FREE(h);                   \
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
    /*SPEX_matrix r,
    SPEX_matrix q,*/
    SPEX_matrix R,    
    SPEX_matrix Q,    
    SPEX_matrix rhos,         // sequence of pivots
    const int64_t j,          // Row of R to compute (col j+1 of Q will also be computed)
    const SPEX_matrix A,      // Matrix to be factored
    int64_t *Prev,
    int64_t *vec,
    const int64_t *leftmost,
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
    int64_t p, pQ, pR, i, top, x,l, prev, iQ;
    int sgn;
    int64_t *xi=NULL, *h=NULL, *col=NULL;

    h = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    // initialize workspace history array
    for (i = 0; i < m; i++)
    {
        h[i] = -1;
    }

    // xi serves as a global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L
    // for the triangular solve.
    xi = (int64_t*) SPEX_malloc(n*sizeof(int64_t)); //TOASK how to free xi becasue we dont use all of it

    col = (int64_t*) SPEX_malloc(n*sizeof(int64_t));

    
    // Obtain the nonzero pattern of the jth row of R (kind of a shitty reach) FIXME
    //printf("j: %ld\n",j);
    l=0;
    for(i=j;i<n;i++)
    {
        for (p=R->p[i];p<R->p[i+1];p++)
        {
            if(R->i[p]==j)
            {
                //printf("i %ld p %ld\n",i,p);
                xi[l]=p;
                col[l]=i;
                l++;
                break;
            }
            else if (R->i[p]>j)
            {
                break;
            }
        }
    }
    
    int64_t estimate = 64 * SPEX_MAX (2, ceil (log2 ((double) n)));
   
    top=0;
    // Compute row k of R
    for (p = top; p < l; p++)
    {
        // Obtain the index of the current nonzero
        x = xi[p];
        i = col[p];//column number
        // R(j,i) = Q(:,j) dot A(:,i)
        //SPEX_MPZ_INIT2(R->x.mpz[x], estimate); //if I do this I get a valgrind leak
        SPEX_CHECK(spex_dot_product(R->x.mpz[x],Q, j, A, i, option)); 
    }
    SPEX_MPZ_SET(rhos->x.mpz[j],R->x.mpz[xi[top]]); //rhos stores the diagonal of R
    SPEX_matrix_check(rhos,option);
    
    // Compute column j+1 of Q using IPGE and history updates (dependent on the j-th column of R)
    for (pQ =Q->p[j+1]; pQ < Q->p[j+2]; pQ++) //Iterate over the nonzeros in col j+1 of Q
    {
        iQ=Q->i[pQ];
        prev=leftmost[iQ];
        
        //if(Prev[pQ]==-1) continue; //if Prev[pQ] is -1, then Q(i,j) is zero and no work must be done
        
            for(pR =R->p[j+1];pR <R->p[j+2];pR++) //Iterate over the nonzeros in col j+1 of R R(:,j+1)
            {
                i=R->i[pR];
                //necesito controlar que i sea igual que el num de columna de Q->x[prev]
                if(i==(j+1)) //row number is the same as column number
                {
                    break;//TODO we know that R has to have a nonzero diagonal, so we can just change the for to pR < R->p[j+2]-1 and avoid lines 117-121
                }
                
                if(prev==pQ) break; 
                if(pQ==10)
                {
                    printf("i %ld, j %ld, h[pQm] %ld\n",i,j+1,h[pQ%m]+1);
                }
                if(i>h[pQ%m]+1) //"an update of Q(i,j)" has been skipped because R(i,l) is zero or Q(i,l) is zero
                {
                    //History update
                    //Q(i,j)=rho^()*Q(i,k)/rho^()
                    // Q[pQ] = x[pQ] * rho[i]
                    printf("hhist i %ld h[pQ] %ld pQ %ld \n",i,h[pQ%m],pQ);
                    SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i-1]);
                    if(i>1 && h[pQ%m]>-1)
                    {
                        SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ%m]]);
                    }
                }

                //check if Q->x.mpz[Prev[pPrev]] is nonzero
                SPEX_MPZ_SGN(&sgn, Q->x.mpz[prev]);
                if(sgn==0) break;                
                    
                SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]); //Q(i,k)==0 this can happen when A(i,k)=0, but Q(i,k) is nonzero
                if(sgn==0)
                {
                    //Q(i,j)=-R(k,j)*Q(i,k)/rho^(k-1)
                    // Q[pQ] = Q[pQ] - R[pR]*Q[k]
                    printf("0 pQ %ld Prev[prev] %ld, prev %ld pR %ld i %ld\n",pQ,Prev[prev],prev,pR,i);

                    SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[prev]); // Q->x.mpz[pQ-m*(j+1-i)
                    if(i>=1)
                    {
                        SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i-1]);
                    }
                }
                else
                {
                    printf("pQ %ld Prev[prev] %ld, prev %ld pR %ld i %ld\n",pQ,Prev[prev],prev,pR,i);
                    //Q(i,j)=(rho^k*Q(i,j)-R(k,j)*Q(i,k))/rho^(k-1)
                    // Q[pQ] = x[pQ] * rho[i]
                    SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i]);
                    // Q[pQ] = Q[pQ] - R[pR]*Q[k]
                    SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[prev]);
                    if(i>=1)
                    {
                        SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i-1]);
                    }
                    
                }
                Prev[vec[iQ]]=pQ;
                vec[iQ]=pQ;
                prev=Prev[prev];

                h[pQ%m]=i;
                
            }
        
        if(h[pQ%m]<j) //to finalize element
        {
            //History update
            //Q(i,j)=rho^()*Q(i,k)/rho^()
            printf("hist h[pQ] %ld j %ld pQ %ld\n",h[pQ%m],j,pQ);
            SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
            if(j>=1 && h[pQ%m]>=0)
            {
                SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ%m]]);
            }
        }

    }
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
