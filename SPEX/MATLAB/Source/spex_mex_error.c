//------------------------------------------------------------------------------
// SPEX/MATLAB/SPEX_mex_error.c: Check error codes for MATLAB
//------------------------------------------------------------------------------

// SPEX: (c) 2022, Chris Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function prints error messages for MATLAB for debugging */

#include "SPEX_mex.h"

void spex_mex_error
(
    SPEX_info status,
    char *message
)
{

    switch (status)
    {
        case SPEX_OK :                  // all is well
            return ;

        case SPEX_OUT_OF_MEMORY :       // out of memory
            SPEX_finalize ( );
            mexErrMsgTxt ("out of memory");

        case SPEX_SINGULAR :            // the input matrix A is singular
            SPEX_finalize ( );
            mexErrMsgTxt ("input matrix is singular");

        case SPEX_INCORRECT_INPUT :     // one or more input arguments are
                                        // incorrect
            SPEX_finalize ( );
            mexErrMsgTxt ("invalid inputs");

//        case SPEX_INCORRECT :         // The solution is incorrect
//            SPEX_finalize ( );
//            mexErrMsgTxt ("result invalid");

        case SPEX_UNSYMMETRIC :         // Matrix is unsymmetric
                                        // cannot use Cholesky
            SPEX_finalize ( );
            mexErrMsgTxt ("input matrix is either pattern or numeric unsymmetric");

        case SPEX_NOTSPD :              // Matrix is not SPD, can't use Cholesky
            SPEX_finalize ( );
            mexErrMsgTxt ("input matrix is symmetric but not SPD");

        case SPEX_PANIC :               // SPEX used without proper
                                        // initialization
            SPEX_finalize ( );
            mexErrMsgTxt ("panic");

        default :
            SPEX_finalize ( );
            mexErrMsgTxt (message);
    }
}

