
#include "SPEX.h"

/* Purpose: Determine why a SPEX_Chol function failed
 */
void SPEX_determine_error
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
    else if (ok == SPEX_UNSYMMETRIC)
    {
        printf("\nInput matrix is unsymmetric, please try left LU\n");
    }
    else if (ok == SPEX_INCORRECT_INPUT)
    {
        printf("\nIncorrect input for a SPEX_Chol Function\n");
    }
}
