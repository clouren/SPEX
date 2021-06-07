//example.h
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

//#include "../../../SuiteSparse_config/SuiteSparse_config.h"
//#include "../../SPEX_Util/Include/SPEX_Util.h"
//#include "../Include/SPEX_Chol.h"
#include "SPEX_Chol.h"

void read(double* input, int size);

int python_backslash
( 
 //output
 double** x,
 //input
 int64_t * col_pointers, 
 int64_t * row_index, 
 double* data,
 double* rhs,
 int n,
 int nz
);

void python_backslash2
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

