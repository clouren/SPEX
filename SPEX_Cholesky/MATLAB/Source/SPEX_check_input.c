//------------------------------------------------------------------------------
// IP_Chol/IP_check_input: Check matlab input
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: This function checks the input of the matlab array. It ensures that
   the input matrix and right hand side vectors are correct format. */
void IP_check_input
(
    const mxArray * input [],     // The matlab input array
    int32_t nargin
)
{
    if (nargin != 2 && nargin != 3)
    {
        mexErrMsgTxt ("incorrect number of inputs") ;
    }

    //--------------------------------------------------------------------------
    // Check matrix A
    //--------------------------------------------------------------------------
    if (mxIsComplex (input[0]))        // Is the matrix real valued?
    {
        mexErrMsgTxt ("Matrix must be real; try backslash instead") ;
    }
    if (!mxIsSparse (input[0]))     // Is the matrix sparse?
    {
        mexErrMsgTxt ("Matrix must be sparse. Try again with A = sparse(A)") ;
    }

    //--------------------------------------------------------------------------
    // Check option
    //--------------------------------------------------------------------------
    if (!mxIsStruct (input[nargin-1]))        // Is third argument the struct?
    {
        mexErrMsgTxt ("Third argument must be the options struct") ;
    }
    if (mxGetNumberOfFields(input[nargin-1]) != 5)
    {
        mexErrMsgTxt("Error! The options struct must have 5 elements. Please "
            "reset it with option = SPEX_get_options;");
    }

    //--------------------------------------------------------------------------
    // Check option
    //--------------------------------------------------------------------------
    if (nargin == 3)
    {
        if (mxIsSparse (input[1]))         // Is b sparse?
        {
            mexErrMsgTxt ("Right hand side vector must be dense. "
                "Try again with b = full(b)") ;
        }
    }
}
