//------------------------------------------------------------------------------
// SPEX_QR/spex_qr_nonzero_structure: Nonzero patten for QR
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

#include "spex_qr_internal.h"


/* Purpose: This function performs a symbolic sparse triangular solve for
 * each column of R
 * It allocates the memory for the R matrix and determines the full nonzero
 * pattern of R
 * It also obtains the nonzero pattern of Q using the column elimination 
 * tree. And fills in the values from A into Q
 *
 * Input arguments of the function:
 *
 * R_handle:    A handle to the R matrix. Null on input.
 *              On output, contains a pointer to the partial RT matrix.
 * 
 * Q_handle:    A handle to the Q matrix. Null on input.
 *              On output, contains a pointer to the partial Q matrix
 *
 * A:           The user's input matrix
 *
 * S:            Symbolic analysis struct for QR factorization.
 *               On input it contains the column permutation
 *               On output it contains the number of nonzeros in R and Q.
 *
 * option:       command options
 */

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
    const SPEX_matrix A,          // Input Matrix
    const SPEX_symbolic_analysis S, // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
    const SPEX_options option
)
{

    // All inputs have been checked by the caller, thus asserts are used here
    // as a reminder of the expected data types
    SPEX_info info;
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);

    int64_t *w, *s, *leftmost;
    int64_t *Qi, *Qp;
    int64_t top, k, len, i, p, n = A->n, m=A->m, m2=m, rnz, qnz, j,h,len2, col,q;
    SPEX_matrix R = NULL, Q=NULL;
    SPEX_matrix QT= NULL, RT=NULL;
    ASSERT(n >= 0);

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    // Allocate R
    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, S->rnz,
        false, false, NULL));

    Qi = (int64_t*) SPEX_malloc((n*m)* sizeof (int64_t));
    Qp = (int64_t*) SPEX_malloc((m+1)* sizeof (int64_t));

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
    // Nonzero pattern of R
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
            
        } 

        //order s
        //qsort(&s[top], n-top, sizeof(int64_t), compare); //TODO check if this is needed, seems like it isn't, but it's staying here in case things start going wrong

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
    R->p[n] = rnz;
    SPEX_CHECK(spex_qr_transpose(&RT,R,NULL));
    //SPEX_matrix_check(R, option);

    
    
    //--------------------------------------------------------------------------
    // Nonzero pattern of QT
    //--------------------------------------------------------------------------
    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */
    
    qnz = 0 ;
    for (k = 0; k < m; k++) //find Q(k,:) pattern
    {    
        Qp [k] = qnz ;  
        top = n ;
        i = leftmost[k];

        for (len = 0 ; i!=-1 && w [i] != k; i = S->parent [i]) /* traverse up to root*/
        {
            s [len++] = i ;
            w [i] = k ;
        }
        while (len > 0) s [--top] = s [--len] ; /* push path on stack */

 
        for (p = top ; p < n ; p++) /* for each i in pattern of Q(:,k) */
        {
            i = s [p] ;                     /* Q(i,k) is nonzero */
            Qi [qnz++] = i ;                  /* Q(i,k) = x(i) */
        }
        
    }
    // Finalize Q->p
    Qp[n] =  qnz;
    SPEX_CHECK(SPEX_matrix_allocate(&QT, SPEX_CSC, SPEX_MPZ, n, m, qnz,
        true, false, NULL));
    
    QT->i = (int64_t*) SPEX_malloc((qnz)* sizeof (int64_t));
    QT->p = (int64_t*) SPEX_malloc((n+1)* sizeof (int64_t));
    memcpy(QT->i, Qi, qnz*sizeof(int64_t));
    memcpy(QT->p, Qp, (n+1)*sizeof(int64_t)); //TODO check m vs n es num cols, pero en matrix m es num cols :/
    QT->p_shallow=false;
    QT->i_shallow=false;

    // Transpose to obtain the nonzero pattern of Q
    SPEX_CHECK(spex_qr_transpose(&Q, QT, NULL));
    Q->nz=qnz; //TODO change
    //printf("n %ld m %ld %ld qnz\n",n,m,qnz );
    //option->print_level = 3;
    //SPEX_matrix_check(Q, option);
    //SPEX_matrix_check(QT, option);
    //--------------------------------------------------------------------------
    // Copy values of A into Q
    //--------------------------------------------------------------------------
    //first column is exactly the same
    for(p=A->p[0];p<A->p[1];p++)
    {
        SPEX_MPZ_SET(Q->x.mpz[p],A->x.mpz[p]);
    }
    //For all other columns the logic is similar to that of the dot product
    for(k=1;k<n;k++) 
    {
        col = S->Q_perm[k];
        p=A->p[col];
        q=Q->p[k];
        while(p < A->p[col+1] && q < Q->p[k+1])
        {
            if(A->i[p] < Q->i[q])
            {
                p++;
            }
            else if(A->i[p] > Q->i[q])
            {
                q++;
            }
            else
            {
                SPEX_MPZ_SET(Q->x.mpz[q],A->x.mpz[p]);
                p++;
                q++;
            }
        }
    }
    
    ///--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------
    (*Q_handle) = Q;
    (*R_handle) = RT; //Return R transpose because of how we store R in factorization


    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
