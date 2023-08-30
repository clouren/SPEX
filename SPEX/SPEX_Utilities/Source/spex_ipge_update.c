//------------------------------------------------------------------------------
// SPEX_Utilities/spex_ipge_update: Performs IPGE update
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_util_internal.h"


SPEX_info spex_ipge_update
(
    SPEX_matrix A,    
    SPEX_matrix B,    
    SPEX_matrix rhos,         // sequence of pivots
    const int64_t a,
    const int64_t b,
    const int64_t i,
    const int64_t j,
    const int64_t k,         
    SPEX_options option
)
{
    
    SPEX_info info;
    
    ASSERT(B->type == SPEX_MPZ);
    ASSERT(B->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(rhos->type == SPEX_MPZ);
    ASSERT(rhos->kind == SPEX_DENSE);
    
    int sgn;
    
    SPEX_MPZ_SGN(&sgn, A->x.mpz[a]);
    if(sgn==0)
    {
        SPEX_MPZ_SUBMUL(A->x.mpz[a], B->x.mpz[b], A->x.mpz[k]);
        if(j>=1)
        {
            SPEX_MPZ_DIVEXACT(A->x.mpz[a], A->x.mpz[a], rhos->x.mpz[i]);
        }
    }
    else
    {
        
        SPEX_MPZ_MUL(A->x.mpz[a], A->x.mpz[a], rhos->x.mpz[j]);
        SPEX_MPZ_SUBMUL(A->x.mpz[a], B->x.mpz[b], A->x.mpz[k]); 
        if(j>=1)
        {
            SPEX_MPZ_DIVEXACT(A->x.mpz[a], A->x.mpz[a], rhos->x.mpz[i]);
        }
    }
    
    return SPEX_OK;
}
