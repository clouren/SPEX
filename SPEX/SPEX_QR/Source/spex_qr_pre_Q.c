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
    int64_t p,q,i, n=A->n, m=A->m, nz=n*m;

    SPEX_CHECK (SPEX_matrix_allocate(&Q, SPEX_CSC, SPEX_MPZ, m, n, m*n, false, false, NULL));
    Q->nz=nz;
    
    //set Q->p
    for(i=0; i<n+1;i++)
    {
        Q->p[i]=i*n;
    }

    for(i=0;i<nz;i++)
    {
        Q->i[i]=i%n;
    }
    
    //this is very inefficient, i don't know how to do thissssss
    for(i=0;i<nz;i++)
    {
        SPEX_MPZ_SET_UI(Q->x.mpz[i],0);
    }

    for(i=0;i<n;i++)
    {
        for(p=A->p[i];p<A->p[i+1];p++)
        {
            for(q=Q->p[i];q<Q->p[i+1];q++)
            {
                //printf("pA %ld iA %ld, pQ %ld iQ %ld\n",p,A->i[p],q,Q->i[p]);
                if(A->i[p]==Q->i[q])
                {
                    SPEX_MPZ_SET(Q->x.mpz[q],A->x.mpz[p]);
                    break;
                }

            }
        }
    }
    
    //SPEX_matrix_check(Q, option);

    (*Q_handle)=Q;
    return SPEX_OK;
}
