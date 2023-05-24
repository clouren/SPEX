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
    SPEX_matrix *Q_handle;
    SPEX matrix A;
    SPEX_options option
)
{
    SPEX_info info;
    
    SPEX_matrix Q;

    SPEX_CHECK (SPEX_matrix_allocate(&Q, SPEX_CSC, SPEX_MPZ, m, n, m*n, false, false, NULL));
    
    //set Q->p
    for(i=0; i<n:i++)
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
        SPEX_MPZ_SETUI(Q->x.mpz[i],0);
    }

    for(i=0;i<n:i++)
    {
        for(p=A->p[i];q<A->p[i+1];p++)
        {
            for(q=Q->p[i];q<Q->p[i+1];p++)
            {
                if(A->i[p]==Q->i[p])
                {
                    SPEX_MPZ_SET(Q->x.mpz[q],A->x.mpz[p]);
                    break;
                }

            }
        }
    }

    (*Q_handle)=Q;
    return SPEX_OK;
}
