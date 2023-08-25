//------------------------------------------------------------------------------
// SPEX_Utilities/spex_history_update: Performs history update
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Lorena Mejia Domenzain, Christopher Lourenco,
// Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_util_internal.h"


SPEX_info spex_history_update
(
    SPEX_matrix A,    
    SPEX_matrix rhos,         // sequence of pivots
    const int64_t a,
    const int64_t i,
    const int64_t j,
    const int64_t k,
    const int64_t thres,         
    SPEX_options option
)
{
    
    SPEX_info info;
    
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(rhos->type == SPEX_MPZ);
    ASSERT(rhos->kind == SPEX_DENSE);
    
    
    //Q(i,j)=rho^()*Q(i,k)/rho^()
    // Q[pQ] = x[pQ] * rho[i]
    SPEX_MPZ_MUL(A->x.mpz[a], A->x.mpz[a], rhos->x.mpz[i]);
    if(j>thres)
    {
        SPEX_MPZ_DIVEXACT(A->x.mpz[a], A->x.mpz[a], rhos->x.mpz[k]);
    }
    
    return SPEX_OK;
}
