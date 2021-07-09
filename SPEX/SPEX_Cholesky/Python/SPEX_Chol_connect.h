//example.h
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

