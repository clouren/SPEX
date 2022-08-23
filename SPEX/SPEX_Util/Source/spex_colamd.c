//------------------------------------------------------------------------------
// SPEX_Util/spex_colamd: Call AMD for matrix ordering
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2022, Chris Lourenco, United States Naval Academy,
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


#include "spex_util_internal.h"

/* Purpose: 
 */
SPEX_info spex_colamd
(
    int64_t **perm_handle,
    int64_t *nnz,
    const SPEX_matrix *A,
    const SPEX_options* option
)
{
    
    SPEX_info info;
    // Check input
    
    int64_t anz; // Number of nonzeros in A
    SPEX_CHECK (SPEX_matrix_nnz(&anz, A, option)) ;
    int64_t i, n = A->n;
    
    int64_t *A2 = NULL , *perm=NULL;
    
    int pr = SPEX_OPTION_PRINT_LEVEL(option);
    
    // Allocate memory for permutation
    perm = (int64_t*)SPEX_malloc( (n+1)*sizeof(int64_t) );
    if (perm == NULL)
    {
        SPEX_FREE_ALL ;
        return (SPEX_OUT_OF_MEMORY) ;
    }
    
    // Declared as per COLAMD documentation
    int64_t Alen = 2*anz + 6 *(n+1) + 6*(n+1) + n;
    A2 = (int64_t*)SPEX_malloc(Alen*sizeof(int64_t));
    if (!A2)
    {
        // out of memory
        SPEX_FREE_ALL ;
        return (SPEX_OUT_OF_MEMORY) ;
    }
    // Initialize S->p as per COLAMD documentation
    for (i = 0; i < n+1; i++)
    {
        perm[i] = A->p[i];
    }
    // Initialize A2 per COLAMD documentation
    for (i = 0; i < anz; i++)
    {
        A2[i] = A->i[i];
    }
    int64_t stats[COLAMD_STATS];
    SuiteSparse_long colamd_result = colamd_l (n, n, Alen,
            (SuiteSparse_long *)A2, (SuiteSparse_long *) perm,
            (double *)NULL, (SuiteSparse_long *) stats);
    if (!colamd_result)
    {
        // COLAMD failed; this "cannot" occur (and is untestable)
        // FIXME make sure this is untestable
        SPEX_FREE_ALL ;
        return (SPEX_PANIC) ;
    }
    // estimate for lnz and unz
    (*nnz) = 10*anz;
    
    // Print stats if desired
    if (pr > 0)
    {
        SPEX_PRINTF ("\n****Ordering Information****\n");
        colamd_l_report ((SuiteSparse_long *) stats);
        SPEX_PRINTF ("\nEstimated L nonzeros: %" PRId64 "\n", nnz);
    }
    SPEX_FREE (A2) ;
    
    (*perm_handle)=perm;
    return SPEX_OK;
}
