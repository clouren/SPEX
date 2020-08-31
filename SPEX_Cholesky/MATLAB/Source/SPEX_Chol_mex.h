//------------------------------------------------------------------------------
// SPEX_Chol/MATLAB/SPEX_Chol_mex.h: include file for MATLAB functions
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_CHOL_MEX_H
#define SPEX_CHOL_MEX_H

#include "SPEX_Chol.h"
#include "matrix.h"

#define SPEX_MEX_OK(method)             \
{                                       \
    status = method ;                   \
    if (status != SPEX_OK)              \
    {                                   \
        spex_chol_mex_error (status, "") ;   \
    }                                   \
}

typedef enum
{
    SPEX_SOLUTION_DOUBLE = 0,   // x as double
    SPEX_SOLUTION_VPA = 1,      // return x as vpa
    SPEX_SOLUTION_CHAR = 2      // x as cell strings
}
SPEX_solution ;

typedef struct spex_mex_options
{
    SPEX_solution solution ;    // how should x be returned to MATLAB
    int32_t digits ;            // # of digits to use for vpa
}
spex_mex_options ;

/* Purpose: A GMP reallocation function
 * This allows GMP to use MATLAB's default realloc function
 */
void* SPEX_gmp_mex_realloc
(
    void* x,    // void* to be reallocated
    size_t a,   // Previous size
    size_t b    // New size
);

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free
 */
void SPEX_gmp_mex_free
(
    void* x,    // void* to be freed
    size_t a    // Size
);

void spex_chol_get_matlab_options
(
    SPEX_options* option,           // Control parameters
    spex_mex_options *mexoptions,   // MATLAB-specific options
    const mxArray* input            // options struct, may be NULL
) ;

// Purpose: This function checks if the array x contains Inf's, NaN's, or
// if its values can be represented as int64_t values.

bool spex_chol_mex_check_for_inf     // true if x can be represented as int64_t
(
    double* x, // The array of numeric values
    mwSize n   // size of array
) ;

/* Purpose: This function reads in the A matrix and right hand side vectors. */
void spex_chol_mex_get_A_and_b
(
    SPEX_matrix **A_handle,     // Internal SPEX Mat stored in CSC
    SPEX_matrix **b_handle,     // mpz matrix used internally
    const mxArray* pargin[],    // The input A matrix and options
    SPEX_options* option
) ;

/* Purpose: Report errors if they arise
 */
void spex_chol_mex_error
(
    SPEX_info status,
    char *message
) ;

#endif

