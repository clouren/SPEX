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
 * col:      Vector containing column numbers
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
/* TOASK
    //Input/Output
    SPEX_matrix R,       // Right triangular matrix
    SPEX_matrix Q,       // Pair-wise orthogonal matrix
    SPEX_matrix rhos,    // sequence of pivots
    int64_t *Qj,         // pointers to elements of the jth column of Q
    int64_t *col,        // column numbers
    int64_t *h,          // History vector
    //Output
    bool *isZeros,       // True if j+1th column of Q is linearly dependent
    //Input
    onst int64_t j,     // Row of R to compute (col j+1 of Q will be finalized)
    const SPEX_matrix A, // Matrix to be factored
    int64_t *Q_perm,     // Column permutation
    SPEX_options option  // Command options
*/

SPEX_info spex_qr_ipgs
(
    SPEX_matrix R,       // Right triangular matrix
    SPEX_matrix Q,       // Pair-wise orthogonal matrix
    SPEX_matrix rhos,    // sequence of pivots
    int64_t *Qj,         // pointers to elements of the jth column of Q
    int64_t *col,        // column numbers
    const int64_t j,     // Row of R to compute (col j+1 of Q will be finalized)
    const SPEX_matrix A, // Matrix to be factored
    int64_t *h,          // History vector
    bool *isZeros,       // True if j+1th column of Q is linearly dependent
    int64_t *Q_perm,     // Column permutation
    SPEX_options option  // Command options
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
    
    *isZeros=true; //start by assuming column of Q is linearly dependent
   
    final = (int64_t*) SPEX_malloc((m)*sizeof(int64_t));
    
    for(k=0;k<m;k++)
    {
        final[k]=-1;
    }

    //--------------------------------------------------------------------------
    // Compute row j of R, store as column
    //--------------------------------------------------------------------------

    for (pR =R->p[j];pR <R->p[j+1];pR++)
    {
        // Obtain the index of the current nonzero
        i = R->i[pR];//column number where j is row number
        // R(j,i) = Q(:,j) dot A(:,i)
        SPEX_CHECK(spex_dot_product(R->x.mpz[pR],Q, j, A, Q_perm[i], option)); 
    }
    //rhos stores the diagonal of R (pivots)
    SPEX_MPZ_SET(rhos->x.mpz[j],R->x.mpz[R->p[j]]);

    //--------------------------------------------------------------------------
    // Update columns j+2 to n of Q (column j+1 is finalized after this)
    //--------------------------------------------------------------------------
    for (pR =R->p[j]; pR < R->p[j+1]; pR++)
    {
        i = R->i[pR];
        if(i<=j+1) continue;//the j+1 column is finalized at the end
    
        //for all elements in column i of Q
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
            //"an update of Q(i,j)" has been skipped because R(i,l) is zero 
            // or Q(i,l) is zero //TODO fix index in comments
            if(j+1>h[pQ]+1)
            {
                SPEX_CHECK(spex_history_update(Q,rhos,pQ,j-1,h[pQ],h[pQ]-1,0,option));
            }
 
            // IPGE update
            SPEX_CHECK(spex_ipge_update(Q,R,rhos,pQ,pR,j-1,j,prev,option));
            
            // Record changes in history vector
            h[pQ]=j+1;
        }
    }
    
    //--------------------------------------------------------------------------
    //for k=j+1
    // IPGE and finalize column j+1 of Q
    //--------------------------------------------------------------------------
    k=j+1;

    // Find the necessary element of R
    for(pR = R->p[j]; pR < R->p[j+1]; pR++)
    {
        i=R->i[pR];
        //printf("i %ld pR %ld \n",i, pR);
        if(i>=k) break; //should happen in the first couple of elements of R[j]
    }

    // For all elements in column k of Q
    for (pQ = Q->p[k]; pQ < Q->p[k+1]; pQ++)
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
                SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], 
                                    rhos->x.mpz[h[pQ]-1]);
            }
            final[iQ]=pQ;
            h[pQ]=k;
        }
        else
        {
            
            if(j+1>h[pQ]+1) //"an u
            {
                //"an update of Q(i,j)" has been skipped because R(i,l) is zero 
                // or Q(i,l) is zero TODO same as before

                SPEX_CHECK(spex_history_update(Q,rhos,pQ,j-1,h[pQ],h[pQ]-1,0,option));
                
            }
            
            if(col[iQ]!=j)
            {
                SPEX_MPZ_MUL(Q->x.mpz[pQ], Q->x.mpz[pQ], rhos->x.mpz[j]);
                if(j>=1)
                {
                    SPEX_MPZ_DIVEXACT(Q->x.mpz[pQ], Q->x.mpz[pQ], 
                                            rhos->x.mpz[j-1]);
                }
                h[pQ]=k;
                continue;
            }
            
            //IPGE update
            SPEX_CHECK(spex_ipge_update(Q,R,rhos,pQ,pR,j-1,j,prev,option));
        }

        // Update history vector
        h[pQ]=k;
        
        // Update the final and col vectors needed for the next iteration
        final[iQ]=pQ;
        col[iQ]=k;
        
        // Checks if the entire column k of Q is zeros for rank revealing QR
        SPEX_MPZ_SGN(&sgn, Q->x.mpz[pQ]);
        if(sgn!=0)
        {
            *isZeros=false;
        }
    }
    
    // Update the final and col vectors needed for the next iteration
    for(i=0;i<m;i++)
    {
        Qj[i]=final[i];
    }
    //SPEX_matrix_check(Q, option); 
    //--------------------------------------------------------------------------
    // Free workspace
    //--------------------------------------------------------------------------
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
