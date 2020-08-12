//------------------------------------------------------------------------------
// SPEX_Left_LU/SPEX_Left_LU_analysis_free: Free memory from symbolic analysis struct
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function frees the memory of the SPEX_LU_analysis struct
 *
 * Input is the SPEX_LU_analysis structure, it is destroyed on function
 * termination.
 */

#include "spex_left_lu_internal.h"

SPEX_info SPEX_Left_LU_analysis_free
(
    SPEX_LU_analysis **S, // Structure to be deleted
    const SPEX_options *option
)
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    if ((S != NULL) && (*S != NULL))
    {
        SPEX_FREE ((*S)->q) ;
        SPEX_FREE (*S) ;
    }

    return (SPEX_OK) ;
}

