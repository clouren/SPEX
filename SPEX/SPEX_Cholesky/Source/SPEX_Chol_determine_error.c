///------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_determine_error: Output error code for SPEX Chol
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Determine what error occured when using IP Chol. This is used solely
 * for demo/debugging purposes.
 */

#include "SPEX_Chol.h"

//TODO Add error code for not symmetric

void SPEX_Chol_determine_error
(
    SPEX_info ok
)
{
    if (ok == SPEX_OUT_OF_MEMORY)
    {
        printf("\nOut of memory\n");
    }
    else if (ok == SPEX_SINGULAR)
    {
        printf("\nInput matrix is singular OR no diagonal pivot. Please ensure input is SPD\n");
    }
    else if (ok == SPEX_INCORRECT_INPUT)
    {
        printf("\nIncorrect input for a IP Chol Function\n");
    }
}
