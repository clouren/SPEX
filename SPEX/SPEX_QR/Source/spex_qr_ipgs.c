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
    SPEX_FREE(col);                   \
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
    /*int64_t *h,         // History vector
    int64_t *xi,        // Nonzero pattern vector of row j of R
    int64_t *c,
    const int64_t *parent,          // Elimination tree*/
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
    int64_t p, pQ, pR, i, top, x,l;
    int sgn;
    int64_t *xi=NULL, *h=NULL, *col=NULL;

    h = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    // initialize workspace history array
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
    }

    // xi serves as a global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L
    // for the triangular solve.
    xi = (int64_t*) SPEX_malloc(n*sizeof(int64_t)); //TOASK how to free xi becasue we dont use all of it

    col = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    // Obtain the nonzero pattern of the jth row of R (kind of a shitty reach)
    l=0;
    for(i=0;i<n;i++)
    {
        for (p=R->p[i];p<R->p[i+1];p++)
        {
            if(R->i[p]==j)
            {
                xi[l]=p;
                col[l++]=i;
                break;
            }
            else if (R->i[p]>j)
            {
                break;
            }
        }
    }

    top=0;
    // Compute row k of R
    for (p = top; p < l; p++)
    {
        // Obtain the index of the current nonzero
        x = xi[p];
        i = col[p];//column number
        // R(j,i) = Q(:,j) dot A(:,i)
        SPEX_MPZ_INIT(R->x.mpz[x]);
        //printf("j %ld x %ld i %ld \n", j,x,i);
        SPEX_CHECK(spex_dot_product(R->x.mpz[x],Q, j, A, i, option)); 
    }
    SPEX_MPZ_SET(rhos->x.mpz[j],R->x.mpz[xi[top]]); //rhos stores the diagonal of R

    // Set Q(:,j)=A(:,j)
    /*for (i = 0; i < m; i++)
    {
        SPEX_mpz_set_ui(q->x.mpz[ i ], 0);
    }
    for (i = A->p[j]; i < A->p[j+1]; i++)
    {
        SPEX_MPZ_SET(q->x.mpz[A->i[i]], A->x.mpz[i]);
    }*/

    // Compute column j+1 of Q using IPGE and history updates (dependent on the j-th column of R)
    for (pQ =Q->p[j+1]; pQ < Q->p[j+2]; pQ++) //if we had a pattern for Q_j this is where it would go
    {
        for(pR =R->p[j+1];pR <R->p[j+2]-1;pR++) //Iterate over the nonzeros in col j+1 of R R(j+1,:)
        {
            i=R->i[pR];
            if(i>h[pQ]+1) //"an update of Q(i,j)" has been skipped 
            {
                //History update
                //Q(i,j)=rho^()*Q(i,k)/rho^()
                // Q[pQ] = x[pQ] * rho[i]
                printf("hist i %ld h[pQ] %ld pQ %ld \n",i,h[pQ],pQ);
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i]);
                if(i>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]]);
                }
            }

            SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]); //Q(i,k)
            if(sgn==0)
            {
                //Q(i,j)=-R(k,j)*Q(i,k)/rho^(k-1)
                // Q[pQ] = Q[pQ] - R[pR]*Q[k]
                SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[pQ-m*(j+1-i)]);
                if(i>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i-1]);
                }
            }
            else
            {
                //Q(i,j)=(rho^k*Q(i,j)-R(k,j)*Q(i,k))/rho^(k-1)
                //DOUBLECHECK, it seems like we don't check for Q(i,j)!=0
                // Q[pQ] = x[pQ] * rho[i]
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i]);
                // Q[pQ] = Q[pQ] - R[pR]*Q[k]
                SPEX_MPZ_SUBMUL(Q->x.mpz[pQ], R->x.mpz[pR], Q->x.mpz[pQ-m*(j+1-i)]);
                if(i>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[i-1]);
                }

                
            }

            h[pQ]=i;
        }

        if(h[pQ]<j) //to finalize element
        {
            //History update
            //Q(i,j)=rho^()*Q(i,k)/rho^()
            printf("hist h[pQ] %ld j %ld pQ %ld\n",h[pQ],j,pQ);
            SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
            if(j>=1 && h[pQ]>=0)
            {
                SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[h[pQ]]);
            }
        }

    }

    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
