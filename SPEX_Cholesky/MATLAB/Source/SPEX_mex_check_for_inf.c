//------------------------------------------------------------------------------
// IP_Chol/IP_mex_check_for_inf: Check A and b for infinity
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"


/* Purpose: This function checks if the input matrix or RHS has numbers too
 * large for double*/
bool IP_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
)
{
 
    
    bool x_is_int64 = true ;

    for (mwSize k = 0; k < n; k++)
    {
        double xk = x [k] ;

        if (mxIsInf (xk))
        {
            mexErrMsgTxt ("A must not have any Inf values") ;
        }

        if (mxIsNaN (xk))
        {
            mexErrMsgTxt ("A must not have any NaN values") ;
        }

        if (x_is_int64)
        {
            if (xk < INT64_MIN || xk > INT64_MAX)
            {
                x_is_int64 = false ;
            }
            else
            {
                int64_t xi = (int64_t) xk ;
                if ((double) xi != xk)
                {
                    x_is_int64 = false ;
                }
            }
        }
    }

    return (x_is_int64) ;
}
