//------------------------------------------------------------------------------
// SPEX_print_options
//------------------------------------------------------------------------------

#include "demos.h"

/* Purpose: This function prints out the user specified/default options.
 * this is primarily intended for debugging
 */

void SPEX_print_options     // display specified/default options to user
(
    SPEX_options option     // options (cannot be NULL) FIXME, why??
)
{

    char *piv, *order;
    if (option->order == SPEX_COLAMD || option->order == SPEX_DEFAULT_ORDERING)
    {
        order = "the COLAMD";
    }
    else if (option->order == SPEX_AMD)
    {
        order = "the AMD";
    }
    else if (option->order == SPEX_NO_ORDERING)
    {
        order = "No";
    }
    else
    {
        order = "(undefined)";
    }

    if (option->pivot == SPEX_SMALLEST)
    {
        piv = "smallest";
    }
    else if (option->pivot == SPEX_DIAGONAL)
    {
        piv = "diagonal";
    }
    else if (option->pivot == SPEX_FIRST_NONZERO)
    {
        piv = "first nonzero";
    }
    else if (option->pivot == SPEX_TOL_SMALLEST)
    {
        piv = "diagonal small tolerance";
    }
    else if (option->pivot == SPEX_TOL_LARGEST)
    {
        piv = "diagonal large tolerance";
    }
    else
    {
        piv = "largest";
    }

    printf("\n\n****COMMAND PARAMETERS****");
    printf("\nUsing %s ordering and selecting the %s pivot", order, piv);
    if (option->pivot == SPEX_TOL_SMALLEST ||
        option->pivot == SPEX_TOL_LARGEST)
    {
        printf("\nTolerance used: %lf\n",option->tol);
    }
}
