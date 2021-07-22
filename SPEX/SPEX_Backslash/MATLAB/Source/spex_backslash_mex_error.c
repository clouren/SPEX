//------------------------------------------------------------------------------
// SPEX_Backslash/MATLAB/spex_backslash_mex_error: Return error messages to matlab
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later
//------------------------------------------------------------------------------

/* Purpose: This function prints error messages for MATLAB for debugging*/

#include "SPEX_Backslash_mex.h"

void spex_backslash_mex_error
(
    SPEX_info status,
    char *message
)
{
    
    switch (status)
    {
        case SPEX_OK :                   // all is well
            return ;

        case SPEX_OUT_OF_MEMORY :        // out of memory
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("out of memory") ;

        case SPEX_SINGULAR :             // the input matrix A is singular
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("input matrix is singular") ;

        case SPEX_INCORRECT_INPUT :      // one or more input arguments are incorrect
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("invalid inputs") ;

        case SPEX_INCORRECT :            // The solution is incorrect
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("result invalid") ;

        case SPEX_PANIC :                // SPEX_LU used without proper initialization
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("panic") ;
        
        default : 
            SPEX_finalize ( ) ;
            mexErrMsgTxt (message) ;
    }
}

