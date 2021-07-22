//------------------------------------------------------------------------------
// SPEX_Cholesky/Python/SPEX_Chol_connect.h: include file to use SPEX_Chol in Python
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#include "SPEX_Chol.h"

void read(double* input, int size);

void python_backslash_char
( 
 //output
 char** sol,                //solution
 //input
 int64_t* col_pointers,     //column pointers of A, an array size is n+1
 int64_t* row_index,        //row indices of A, of size nzmax.
 double* data,              //values of A
 double* rhs,               //values of b
 int n,                     // # of rows/columns of A
 int nz                     // max # entries for CSC matrix A
);

void python_backslash_double
( 
 //output
 double* sol,                //solution
 //input
 int64_t* col_pointers,     //column pointers of A, an array size is n+1
 int64_t* row_index,        //row indices of A, of size nzmax.
 double* data,              //values of A
 double* rhs,               //values of b
 int n,                     // # of rows/columns of A
 int nz                     // max # entries for CSC matrix A
);

