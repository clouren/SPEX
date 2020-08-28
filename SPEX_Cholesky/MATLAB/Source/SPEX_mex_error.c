//------------------------------------------------------------------------------
// IP_Chol/IP_mex_error: Return error to matlab
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: This function prints error messages for MATLAB */
void IP_mex_error
(
    SPEX_info status
)
{
    if (status == SPEX_OUT_OF_MEMORY)
    {
        mexErrMsgTxt("Error, Out of memory");
    }
    else if (status == SPEX_SINGULAR)
    {
        mexErrMsgTxt("Error, Input matrix is singular");
    }
    else if (status == SPEX_INCORRECT_INPUT)
    {
        mexErrMsgTxt("Error, Input is incorrect");
    }
}
