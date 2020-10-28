//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR.h: user #include file for SPEX_QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020, Chris Lourenco, United States Naval Academy. 
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_QR_H
#define SPEX_QR_H


// This software performs an exact integer-preserving QR factorization
// WARNING: This code is experimental and developmental, please do not use it.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Chris Lourenco (chrisjlourenco@gmail.com, lourenco@usna.edu)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX_QR is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Christopher Lourenco
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX QR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    To be done


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX QR------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SPEX_Util.h"

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_QR_VERSION "0.0.1"
#define SPEX_CHOL_VERSION_MAJOR 0
#define SPEX_CHOL_VERSION_MINOR 0
#define SPEX_CHOL_VERSION_SUB   1


/* Compute the dot product of two integer vectors x,y and return in z */
void SPEX_dot
(
    SPEX_matrix* x,
    SPEX_matrix* y,
    mpz_t z
)
{
    // Check inputs, x and y must both be column vectors of identical size and
    // stored as dense
    ASSERT( x->n == y->n);
    ASSERT( x->m == y->m);
    ASSERT( x->type == SPEX_MPZ);
    ASSERT( y->type == SPEX_MPZ);
    ASSERT( x->kind == SPEX_DENSE);
    ASSERT( y->kind == SPEX_DENSE);
    
    int64_t k;
    
    // Set z = 0
    SPEX_mpz_set_ui(z, 0);
    for (k = 0; k < x->m; k++)
    {
        SPEX_mpz_addmul(z, x->x.mpz[k], y->x.mpz[k]);
    }
}

// Much more to come

#endif
