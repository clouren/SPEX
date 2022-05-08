//------------------------------------------------------------------------------
// SPEX_Util/SPEX_determine_symmetry: Determine if given matrix is 
//                                    *numerically* (thus pattern-wise) symmetric
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: Determine if the input A is *numerically* (thus pattern-wise) symmetric.
 * Since SPEX is an exact framework, it doesn't make sense to check only pattern symmetry.
 * 
 * If the matrix is determined to be symmetric, SPEX_OK is returned; otherwise,
 * SPEX_UNSYMMETRIC is returned.
 */


// TODO: propagate the new calling style throughout

#define SPEX_FREE_ALL               \
{                                   \
    SPEX_matrix_free(&T,NULL);      \
    SPEX_matrix_free(&R,NULL);      \
}

#include "spex_util_internal.h"

SPEX_info SPEX_determine_symmetry
(
    SPEX_matrix* A,            // Input matrix to be checked for symmetry
    const SPEX_options* option // Command options
)
{    
    SPEX_info info;
    // TODO INSERT UNINITIALIZED CRAP
  
    // Only used index
    int64_t j;

    // Declare matrices T and R. T = A' and R = T' = A''
    SPEX_matrix *T = NULL, *R = NULL ;
    // T = A'
    SPEX_CHECK( SPEX_transpose(&T, A, option) );

    // Check if the number of nonzeros in the columns
    // of A are equal to the number of nonzeros in 
    // the rows of A. This is a quick check to 
    // ensure the matrix is candidate to be symmetric.
    // Moreover, this check is important becuase 
    // otherwise the ensuing block could seg-fault :(
    for (j = 0; j <= A->n; j++)
    {
        if (T->p[j] != A->p[j])
        {
            // nnz( A(:,k)) != nnz( A(k,:))
            SPEX_FREE_ALL ;
            return SPEX_UNSYMMETRIC;
        }
      
    }

    // Set R = T'
    SPEX_CHECK( SPEX_transpose(&R, T, option) );
    
    // Check whether A[i][j] = A[j][i] in both pattern and numerics
    for (j = 0; j < R->p[R->n]; j++)
    {
        // Check pattern
        if (T->i[j] != R->i[j])
        {
            // Not pattern symmetrc
            SPEX_FREE_ALL ;
            return SPEX_UNSYMMETRIC;
        }
        // Check numerics
        int r;
        SPEX_CHECK( SPEX_mpz_cmp(&r, R->x.mpz[j], T->x.mpz[j]) );
        if (r != 0)
        {
            // Not numeric symmetric
            SPEX_FREE_ALL ;
            return SPEX_UNSYMMETRIC;
        }
    }
    
    // Free memory and return OK meaning the matrix is symmetric
    SPEX_FREE_ALL ;
    return SPEX_OK;
}

