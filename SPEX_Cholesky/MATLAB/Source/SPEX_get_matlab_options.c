//------------------------------------------------------------------------------
// IP_Chol/IP_get_matlab_options: Set options for matlab cholesky
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: This function reads in the necessary information from the options
   struct */
void IP_get_matlab_options
(
    SPEX_options* option,  // Control parameters
    const mxArray* input   // The input options from MATLAB interface
)
{
    mxArray* tmp;

    // Get the column ordering
    tmp = mxGetField(input, 0, "column");
    if (tmp == NULL)
    {
        mexErrMsgTxt("Error for column ordering");
    }
    int32_t order = (int32_t) mxGetScalar(tmp);

    // Get the pivoting scheme
    tmp = mxGetField(input, 0, "pivot");
    if (tmp == NULL)
    {
        mexErrMsgTxt("Error at getting pivot");
    }
    int32_t piv = (int32_t) mxGetScalar(tmp);

    // Tolerance if some form of tolerance partial pivoting is used
    if (piv == 3 || piv == 4)
    {
        tmp = mxGetField(input, 0, "tol");
        if (tmp == NULL)
        {
            mexErrMsgTxt("Error at getting the tolerance parameter");
        }
        option->tol = mxGetScalar(tmp);
    }

    //--------------------------------------------------------------------------
    // Verify that the parameters are correct
    //--------------------------------------------------------------------------
    if (order <= 2 && order >= 0) {option->order = (SPEX_col_order) order;}
    if (piv <= 5 && piv >= 0) {option->pivot = (SPEX_pivot) piv;}
    if (option->tol > 1 || option->tol <= 0) {option->tol = 0.1;}
}
