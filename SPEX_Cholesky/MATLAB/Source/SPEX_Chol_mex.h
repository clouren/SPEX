//------------------------------------------------------------------------------
// IP-Chol/MATLAB/IP-Chol_mex.h: Include file for IP-Chol
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#ifndef IP_mex 
#define IP_mex

#include "../../Include/IP-Chol.h"
# include "mex.h"
#include "matrix.h"


#define IP_MEX_OK(method)         \
{                                   \
    status = method;                \
    IP_mex_error(status);         \
}


/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */
void* IP_gmp_mex_realloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
);

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free 
 */
void IP_gmp_mex_free 
(
    void* x,    // void* to be freed
    size_t a    // Size
);

void IP_check_input
(
    const mxArray * input [],    // The matlab input array
    int32_t nargin
);

void IP_get_matlab_options
(
    SLIP_options* option,  // Control parameters
    const mxArray* input   // The input options from MATLAB interface
);


/* Purpose: Check the input x array for numbers too large for 
 * double precision.
 */
bool IP_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
);

/* Purpose: This function reads in the A matrix and right hand side vectors. */
void IP_mex_get_A_and_b
(
    SLIP_matrix **A_handle,  // Internal SLIP Mat stored in CSC
    SLIP_matrix **b_handle,   // mpz matrix used internally
    const mxArray* pargin[], // The input A matrix and options
    int nargin,               // Number of input to the mexFunction
    SLIP_options* option
);


/* Purpose: Report errors if they arise
 */
void IP_mex_error
(
    SLIP_info status
) ;


#endif
