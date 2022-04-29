//------------------------------------------------------------------------------
// SPEX_Left_LU/MATLAB/spex_left_lu_mex_error: Return error messages to matlab
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.
//------------------------------------------------------------------------------

/* Purpose: This function prints error messages for MATLAB for debugging*/

#include "SPEX_Chol_mex.h"

void spex_chol_mex_error
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
            
        case SPEX_UNSYMMETRIc:          // the input matrix A is unsymmetric
            SPEX_finalize();
            mexErrMsgTxt("input matrix is unsymmetric");
            
        case SPEX_NOTSPD :             // the input matrix A is not spd
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("input matrix is not spd") ;


        case SPEX_PANIC :                // SPEX_LU used without proper initialization
            SPEX_finalize ( ) ;
            mexErrMsgTxt ("panic") ;
        
        default : 
            SPEX_finalize ( ) ;
            mexErrMsgTxt (message) ;
    }
}

