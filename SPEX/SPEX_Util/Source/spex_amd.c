//------------------------------------------------------------------------------
// SPEX_Util/spex_amd: Call AMD for matrix ordering
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2022, Chris Lourenco, United States Naval Academy,
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_ALL                           \
{                                               \
    SPEX_FREE_WORKSPACE ;                       \
    SPEX_free(&perm_handle, option);            \
}
//TODO check this free is correct with valgrind

#include "spex_util_internal.h"

/* Purpose: 
 */
SPEX_info spex_amd
(
    int64_t **perm_handle,
    int64_t *nnz,
    const SPEX_matrix *A,
    const SPEX_options* option
)
{
    
    int pr = SPEX_OPTION_PRINT_LEVEL(option);
    int64_t n = A->n;
    int64_t *perm=NULL;
    
    perm = (int64_t*)SPEX_malloc( (n+1)*sizeof(int64_t) );
    if (perm == NULL)
    {
        SPEX_FREE_ALL ;
        return (SPEX_OUT_OF_MEMORY) ;
    }
    
    double Control[AMD_CONTROL];           // Declare AMD control
    amd_l_defaults(Control);              // Set AMD defaults
    double Info [AMD_INFO];
    // Perform AMD
    SuiteSparse_long amd_result = amd_l_order(n,
                (SuiteSparse_long *)A->p, (SuiteSparse_long *)A->i,
                (SuiteSparse_long *)perm, Control, Info);
    if (pr > 0)   // Output AMD info if desired
    {
        SPEX_PRINTF("\n****Ordering Information****\n");
        amd_l_control(Control);
        amd_l_info(Info);
    }
    if (!(amd_result == AMD_OK || amd_result == AMD_OK_BUT_JUMBLED))
    {
        // AMD failed: either out of memory, or bad input
        SPEX_FREE_ALL ;
        if (amd_result == AMD_OUT_OF_MEMORY)
        {
            // AMD ran out of memory
            return (SPEX_OUT_OF_MEMORY) ;
        }
        // input matrix is invalid
        return (SPEX_INCORRECT_INPUT) ;
    }
    (*nnz) = Info[AMD_LNZ];  // Exact number of nonzeros for Cholesky
    
    (*perm_handle)=perm;
    return SPEX_OK;
}
