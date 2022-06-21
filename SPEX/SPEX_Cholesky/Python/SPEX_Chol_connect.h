//------------------------------------------------------------------------------
// SPEX_Cholesky/Python/SPEX_Chol_connect.h: include file to use SPEX_Chol in Python
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "SPEX.h"


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
