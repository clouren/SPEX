//------------------------------------------------------------------------------
// SPEX_Cholesky/Python/SPEX_Chol_connect.h: include file to use SPEX_Chol in Python
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------
/*
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>*/

#include "SPEX.h"


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
 int nz,                     // max # entries for CSC matrix A
 int ordering               // type of ordering: 0-none, 1-colamd, 2-amd
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
 int nz,                     // max # entries for CSC matrix A
 int ordering               // type of ordering: 0-none, 1-colamd, 2-amd
);

/*
//------------------------------------------------------------------------------
// SPEX_Cholesky/Python/SPEX_Chol_connect.h: include file to use SPEX_Chol in Python
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// TODO Delete me?
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#include "SPEX_Chol.h"
*/

SPEX_info SPEX_python_backslash
( 
     //output
     void** sol_void, //solution
     //input
     int64_t* Ap,     // column pointers of A, an array size is n+1 
     int64_t* Ai,     // row indices of A, of size nzmax.
     double* Ax,      // values of A
     double* bx,      // values of b
     int m,           // Number of rows of A
     int n,           // Number of columns of A
     int nz,          // Number of nonzeros in A
     int ordering,    // type of ordering: 0-none, 1-colamd, 2-amd
     bool charOut    // True if char** output, false if double
);

SPEX_info SPEX_python_backslashVoid
( 
//output
 char** sol_char,      // solution (null if double output)
 double* sol_doub,     // solution (null if char output)
 void** sol_void,
 //input
 bool charOut,    // True if char** output, false if double
 int64_t* Ap,     // column pointers of A, an array size is n+1 
 int64_t* Ai,     // row indices of A, of size nzmax.
 double* Ax,      // values of A
 double* bx,      // values of b
 int m,           // Number of rows of A
 int n,           // Number of columns of A
 int nz,           // Number of nonzeros in A
 int ordering     // type of ordering: 0-none, 1-colamd, 2-amd
);
/**/
