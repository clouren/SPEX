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


SPEX_info spex_qr_pre_Q
(
    SPEX_matrix *Q_handle,
    SPEX_matrix A,
    SPEX_options option
)
{
    SPEX_info info;
   
    SPEX_matrix Q=NULL;
    int64_t p,q,i,j, n=A->n, m=A->m, nz=n*m;

    SPEX_CHECK (SPEX_matrix_allocate(&Q, SPEX_CSC, SPEX_MPZ, m, n, m*n, false, true, NULL));
    Q->nz=nz;
    
    //set Q->p
    /*for(i=0; i<n+1;i++)
    {
        Q->p[i]=i*m;
    }

    for(i=0;i<nz;i++)
    {
        Q->i[i]=i%m;
    }
    int64_t estimate = 64 * SPEX_MAX (2, ceil (log2 ((double) n)));
    
    //this is very inefficient, i don't know how to do thissssss
    for(i=0;i<nz;i++)
    {
        SPEX_MPZ_INIT2(Q->x.mpz[i], estimate);
    }
    
   for(j=0;j<n;j++) //this works because Q is "dense"
    {
        for(p=A->p[j];p<A->p[j+1];p++)
        {
            i=A->i[p];
            if(i==Q->i[j*n+i])
            {
                SPEX_MPZ_SET(Q->x.mpz[j*n+i],A->x.mpz[p]);
            }

        }
    }*/
    
    //compute tree of A
    
    
    //get the rows of Q using the tree of A
   
    
    //SPEX_matrix_check(Q, option);

    (*Q_handle)=Q;
    return SPEX_OK;
}
