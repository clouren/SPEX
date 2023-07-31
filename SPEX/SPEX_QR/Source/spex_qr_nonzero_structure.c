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
}

# define SPEX_FREE_ALL               \
{                                    \
    SPEX_FREE_WORKSPACE              \
    SPEX_matrix_free(&L, NULL);      \
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
    const SPEX_symbolic_analysis S  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
)
{

    // All inputs have been checked by the caller, thus asserts are used here
    // as a reminder of the expected data types
    SPEX_info info;
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);

    //int64_t  top, k, j, jnew, n = A->n, p = 0;
    int64_t *w, *s, *leftmost;
    int64_t top, k, len, i, p, n = A->n, m2=n, m=A->m, rnz, qnz, j,h,len2;
    //int64_t *c = NULL;
    SPEX_matrix R = NULL, Q=NULL;
    SPEX_matrix QT= NULL, L=NULL;
    ASSERT(n >= 0);

    //--------------------------------------------------------------------------
    // Declare memory for L and c
    //--------------------------------------------------------------------------

    // Allocate R
    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, S->unz,
        false, false, NULL));
    SPEX_CHECK(SPEX_matrix_allocate(&L, SPEX_CSC, SPEX_MPZ, n, n, S->unz,
        false, false, NULL));

    SPEX_CHECK(SPEX_matrix_allocate(&QT, SPEX_CSC, SPEX_MPZ, m, n, m*n,
        false, false, NULL));

 
    // Allocate c
    //c = (int64_t*) SPEX_malloc(n* sizeof (int64_t));

    w = (int64_t*) SPEX_malloc((n+m2)* sizeof (int64_t));
    leftmost = (int64_t*) SPEX_malloc(m* sizeof (int64_t));
    s = w + n ;
    if (!w)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    
    

    R->i[0] = 0;
    //c[0]++;


    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */

    for (i = 0 ; i < m ; i++) leftmost [i] = -1 ;
    for (k = n-1 ; k >= 0 ; k--)
    {
        for (p = A->p [k] ; p < A->p [k+1] ; p++)
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

        for (p = A->p [k] ; p < A->p [k+1] ; p++)   /* find R(:,k) pattern */
        {
            i = leftmost [A->i [p]] ;         /* i = min(find(A(i,q))) */
            for (len = 0 ; w [i] != k ; i = S->parent [i]) /* traverse up to k */
            {
                s [len++] = i ;
                w [i] = k ;
            }

            while (len > 0) s [--top] = s [--len] ; /* push path on stack */
        }
    
        //order s
        bubbleSort(s,top,n);
        

        for (p = top ; p < n ; p++) /* for each i in pattern of R(:,k) */
        {
            i = s [p] ;                     /* R(i,k) is nonzero */
            R->i [rnz++] = i ;                  /* R(i,k) = x(i) */
        }
        R->i [rnz++] = k ;                     /* R(k,k) */
    }
    // Finalize R->p
    R->p[n] = S->unz = rnz;
    SPEX_CHECK(spex_qr_transpose(&L,R,NULL));
    (*R_handle) = L;//this is a messs, FIXME
    /*

    SPEX_options option = NULL;
    SPEX_create_default_options(&option);
    option->print_level = 3;
    SPEX_matrix_check(R, option);*/

    // Q
    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */

    
    qnz = 0 ;
    for (k = 0; k < n; k++) //find Q(k,:) pattern
    {
        QT->p [k] = qnz ;      

        top = n ;
        p=A->p[k];

            i = leftmost[k];
            //printf("k %ld  left %ld\n",k,i);
            for (len = 0 ; i!=-1 && w [i] != k; i = S->parent [i]) /* traverse up to root*/
            {
                //printf("k %ld p %ld i: %ld\n",k,p,i);
                s [len++] = i ;
                w [i] = k ;
            }
            while (len > 0) s [--top] = s [--len] ; /* push path on stack */

 
        for (p = top ; p < n ; p++) /* for each i in pattern of Q(:,k) */
        {
            i = s [p] ;                     /* Q(i,k) is nonzero */
            QT->i [qnz++] = i ;                  /* Q(i,k) = x(i) */
        }
        
    }
    // Finalize Q->p
    QT->p[n] = S->lnz = qnz;
    SPEX_CHECK(spex_qr_transpose(&Q, QT, NULL));
    Q->nz=qnz;
    //copy A into Q
    //first column is exactly the same
    for(p=A->p[0];p<A->p[1];p++)
    {
        SPEX_MPZ_SET(Q->x.mpz[p],A->x.mpz[p]);
    }
    
    for(k=1;k<n;k++)
    {
        h=0; //h makes it so that we don't start at the begining of the column of Q every single time, but we start where we left off
        for(p=A->p[k];p<A->p[k+1];p++)
        {
            i=A->i[p];
            for(j=Q->p[k]+h;j<Q->p[k+1];j++)
            {
                if(i==Q->i[j])
                {
                    h=j-Q->p[k];
                    SPEX_MPZ_SET(Q->x.mpz[j],A->x.mpz[p]);
                    continue;
                }
                
            }
        }
        
    }


    (*Q_handle) = Q;
    
    /*SPEX_options option = NULL;
    SPEX_create_default_options(&option);
    option->print_level = 3;
    SPEX_matrix_check(A, option);
    SPEX_matrix_check(Q, option);*/


    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
