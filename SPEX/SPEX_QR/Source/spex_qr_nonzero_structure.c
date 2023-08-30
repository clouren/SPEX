//------------------------------------------------------------------------------
// SPEX_QR/spex_qr_pre_factor: Symbolic left-looking Cholesky for QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Lorena Mejia Domenzain, Christopher Lourenco,
// Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_FREE(w);                   \
    SPEX_FREE(leftmost);            \
    SPEX_matrix_free(&QT,NULL);     \
    SPEX_matrix_free(&R, NULL);      \
    SPEX_FREE(Qi);                   \
    SPEX_FREE(Qp);                   \
}

# define SPEX_FREE_ALL               \
{                                    \
    SPEX_FREE_WORKSPACE              \
    SPEX_matrix_free(&RT, NULL);      \
    SPEX_matrix_free(&Q,NULL);       \
}

#include "spex_cholesky_internal.h"


/* Purpose: This function performs a symbolic sparse triangular solve for
 * each column of R
 * It allocates the memory for the R matrix and determines the full nonzero
 * pattern of R
 * It also obtains the nonzero pattern of Q using the column elimination 
 * tree.
 *
 * Importantly, this function assumes that A has already been permuted.
 *
 * Input arguments of the function:
 *
 * R_handle:    A handle to the R matrix. Null on input.
 *              On output, contains a pointer to the partial R matrix.
 * 
 * Q_handle:    A handle to the Q matrix. Null on input.
 *              On output, contains a pointer to the partial Q matrix
 *
 * A:           The user's permuted input matrix
 *
 * S:            Symbolic analysis struct for QR factorization.
 *               On input it contains information that is not used in this
 *               function such as the row/column permutation
 *               On output it contains the number of nonzeros in R.
 */
void swap(int64_t* xp, int64_t* yp)
{
    int64_t temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// An optimized version of Bubble Sort
void bubbleSort(int64_t *arr, int64_t top,int64_t n)
{
    int i, j;
    bool swapped;
    
    for (i = 0; i < n - 1; i++) 
    {
        swapped = false;
        for (j = top; j < n - i - 1; j++) 
        {
            if (arr[j] > arr[j + 1]) 
            {
                swap(&arr[j], &arr[j + 1]);
                swapped = true;
            }
        }
        // If no two elements were swapped by inner loop,
        // then break
        if (swapped == false) break;
        
    }
    
}

// Sorting function
static inline int compare (const void * a, const void * b)
{
    return ( *(int64_t*)a - *(int64_t*)b );
}


SPEX_info spex_qr_nonzero_structure
(
    // Output
    SPEX_matrix *R_handle,        // On output: partial R matrix
                                  // On input: undefined
    SPEX_matrix *Q_handle,        // On output: partial R matrix
                                  // On input: undefined
    // Input
    //int64_t *xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix A,          // Input Matrix
    const SPEX_symbolic_analysis S,  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
    SPEX_options option
)
{

    // All inputs have been checked by the caller, thus asserts are used here
    // as a reminder of the expected data types
    SPEX_info info;
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);

    int64_t *w, *s, *leftmost;
    int64_t *Qi, *Qp;
    int64_t top, k, len, i, p, n = A->n, m=A->m, m2=m, rnz, qnz, j,h,len2, col;
    SPEX_matrix R = NULL, Q=NULL;
    SPEX_matrix QT= NULL, RT=NULL;
    ASSERT(n >= 0);

    //--------------------------------------------------------------------------
    // Declare memory for R
    //--------------------------------------------------------------------------

    // Allocate R
    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, S->rnz,
        false, false, NULL));

    Qi = (int64_t*) SPEX_malloc((n*m)* sizeof (int64_t));
    Qp = (int64_t*) SPEX_malloc((n*m)* sizeof (int64_t));

    w = (int64_t*) SPEX_malloc((n+m2)* sizeof (int64_t));
    leftmost = (int64_t*) SPEX_malloc(m* sizeof (int64_t));
    s = w + n ;
    if (!w)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    
    

    R->i[0] = 0;


    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */

    for (i = 0 ; i < m ; i++) leftmost [i] = -1 ;
    for (k = n-1 ; k >= 0 ; k--)
    {
        col = S->Q_perm[k];
        for (p = A->p [col] ; p < A->p [col+1] ; p++)
        {
            leftmost [A->i [p]] = k ;         /* leftmost[i] = min(find(A(i,:)))*/
        }
    }

    //--------------------------------------------------------------------------
    // Iterations 1:n-1
    //--------------------------------------------------------------------------
    rnz = 0 ;
    for (k = 0; k < n; k++)
    {
        R->p [k] = rnz ;      
        w [k] = k ;  
        top = n ;
        col = S->Q_perm[k];

        for (p = A->p [col] ; p < A->p [col+1] ; p++)   /* find R(:,k) pattern */
        {
            i = leftmost [A->i [p]] ;         /* i = min(find(A(i,q))) */
            for (len = 0 ; w [i] != k ; i = S->parent [i]) /* traverse up to k */
            {
                s [len++] = i ;
                w [i] = k ;
            }
            while (len > 0) s [--top] = s [--len] ; /* push path on stack */
            
            i=A->i [p];
            //i=pinv[A->i [p]];
            if (i > k && w [i] < k)         /* pattern of V(:,k) = x (k+1:m) */
            {
                w [i] = k ;
            }
        } 

        //order s
        bubbleSort(s,top,n);//FIXME qsort(&xi[top], n-top, sizeof(int64_t), compare);     

        for (p = top ; p < n ; p++) /* for each i in pattern of R(:,k) */
        {
            i = s [p] ;                     /* R(i,k) is nonzero */
            R->i [rnz] = i ;                  /* R(i,k) = x(i) */
            rnz++;
        }
        R->i [rnz] = k ;                     /* R(k,k) */
        rnz++;
    }
    // Finalize R->p
    R->p[n] = S->rnz = rnz;
    SPEX_CHECK(spex_qr_transpose(&RT,R,NULL));
    (*R_handle) = RT;


    // Q
    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */
    
    qnz = 0 ;
    for (k = 0; k < m; k++) //find Q(k,:) pattern
    {    
        Qp [k] = qnz ;  
        top = n ;
        col = S->Q_perm[k];

        //i = leftmost[col];
        i = leftmost[k];
        //printf("k %ld col %ld  left %ld\n",k,col,i);
        for (len = 0 ; i!=-1 && w [i] != k; i = S->parent [i]) /* traverse up to root*/
        {
            //printf("k %ld p %ld i: %ld\n",k,p,i);
            s [len++] = i ;
            w [i] = k ;
        }
        while (len > 0) s [--top] = s [--len] ; /* push path on stack */

 
        for (p = top ; p < n ; p++) /* for each i in pattern of Q(:,k) */
        {
            //printf("p %ld i %ld\n", p, s[p]);
            i = s [p] ;                     /* Q(i,k) is nonzero */
            Qi [qnz++] = i ;                  /* Q(i,k) = x(i) */
        }
        
    }
    // Finalize Q->p
    Qp[n] = S->qnz = qnz;
    SPEX_CHECK(SPEX_matrix_allocate(&QT, SPEX_CSC, SPEX_MPZ, n, m, S->qnz,
        true, false, NULL));
    /*bool ok ;
    Qi = (int64_t *) SPEX_realloc (qnz, n*m, sizeof (int64_t), Qi, &ok);
    if (!ok)    {return SPEX_OUT_OF_MEMORY;}
    Qp = (int64_t *) SPEX_realloc (qnz, n*m, sizeof (int64_t), Qp, &ok);
    if (!ok)    {return SPEX_OUT_OF_MEMORY;}*/
    QT->i = (int64_t*) SPEX_malloc((qnz)* sizeof (int64_t));
    QT->p = (int64_t*) SPEX_malloc((qnz)* sizeof (int64_t));
    memcpy(QT->i, Qi, qnz*sizeof(int64_t));
    memcpy(QT->p, Qp, qnz*sizeof(int64_t));
    QT->p_shallow=false;
    QT->i_shallow=false;
    SPEX_CHECK(spex_qr_transpose(&Q, QT, NULL));
    Q->nz=qnz;

    //copy A into Q
    //first column is exactly the same
    for(p=A->p[0];p<A->p[1];p++)
    {
        SPEX_MPZ_SET(Q->x.mpz[p],A->x.mpz[p]);
    }
    
    for(k=1;k<n;k++) //FIXME there has to be a better more sparse way of doing this, maybe look at dotprod
    {
        h=0; //h makes it so that we don't start at the begining of the column of Q every single time, but we start where we left off
        col = S->Q_perm[k];
        for(p=A->p[col];p<A->p[col+1];p++)
        {
            i=A->i[p];
            for(j=Q->p[col]+h;j<Q->p[col+1];j++)
            {
                if(i==Q->i[j])
                {
                    h=j-Q->p[col];
                    SPEX_MPZ_SET(Q->x.mpz[j],A->x.mpz[p]);
                    continue;
                }
                
            }
        }
        
    }


    (*Q_handle) = Q;

    for (i = 0; i < n; i++)
    {
        printf("%ld %ld: ", i, Q->p[i]);
        for(p=Q->p[i];p<Q->p[i+1];p++)
        {
            printf("%ld %ld, ", p, Q->i[p]);
        }
        printf("\n");
    }

    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
