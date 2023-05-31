//------------------------------------------------------------------------------
// SPEX_QR/spex_qr_permute_A: Permutation of matrix A
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Lorena Mejia Domenzain, Christopher Lourenco,
// Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

//TOASK this feels a little silly because it is basically the same as cholesky_permute_A except for Q_perm instead of P_perm

#include "spex_qr_internal.h"

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL { SPEX_matrix_free (&PAQ, NULL); }

/* Purpose: Given the column permutation Q and an inverse row permutation Pinv 
 * stored in S, permute the matrix A and return PAQ'
 * Input arguments:
 *
 * PAQ_handle:   The user's permuted input matrix.
 *
 * A:            The user's input matrix
 *
 * S:            Symbolic analysis struct for Cholesky factorization.
 *               Contains row/column permutation of A
 */

SPEX_info spex_qr_permute_A
(
    //Output
    SPEX_matrix* PAQ_handle,   // On input: undefined
                               // On output: contains the permuted matrix
    //Input
    const SPEX_matrix A,       // Input matrix
    const bool numeric,        // True if user wants to permute pattern and
                               // numbers, false if only pattern
    const SPEX_symbolic_analysis S,  // Symbolic analysis struct that contains
                                // row/column permutations
    const SPEX_options option   // Command options (Default if NULL)
)
{

    SPEX_info info;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    ASSERT(A != NULL);
    ASSERT(S != NULL);
    ASSERT(PAQ_handle != NULL);
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);

    // Create indices and pinv, the inverse row permutation
    int64_t j, k, t, nz = 0, n = A->n, m=A->m;
    (*PAQ_handle) = NULL ;
    //int64_t *pinv = NULL;

    // Allocate memory for PAQ which is a permuted copy of A
    SPEX_matrix PAQ = NULL ;
    SPEX_CHECK(SPEX_matrix_allocate(&PAQ, SPEX_CSC, SPEX_MPZ, m, n, A->p[n],
        false, true, NULL));

    
    if(numeric)
    {

        //----------------------------------------------------------------------
        // construct PAQ with numerical values
        //----------------------------------------------------------------------

        // Set PAQ scale
        SPEX_MPQ_SET(PAQ->scale, A->scale);

        // Populate the entries in PAQ
        for (k = 0; k < n; k++)
        {
            // Set the number of nonzeros in the kth column of PAQ
            PAQ->p[k] = nz;
            // Column k of PAQ is equal to column S->Q_perm[k] of A. j is the
            // starting point for nonzeros and indices for column S->Q_perm[k]
            // of A
            j = S->Q_perm[k];
            // Iterate across the nonzeros in column S->P_perm[k]
            for (t = A->p[j]; t < A->p[j+1]; t++)
            {
                // Set the nonzero value and location of the entries in column
                // k of PAQ
                SPEX_MPZ_SET(PAQ->x.mpz[nz], A->x.mpz[t]);
                // Row i of this nonzero is equal to A->i[t]
                PAQ->i[nz] = A->i[t];
                // Move to the next nonzero element of PAQ
                nz++;
            }
        }
    }
    else
    {

        //----------------------------------------------------------------------
        // construct PAQ with just its pattern, not the values
        //----------------------------------------------------------------------

        SPEX_FREE (PAQ->x.mpz);
        ASSERT (PAQ->x.mpz == NULL);
        PAQ->x_shallow = true ;

        // Populate the entries in PAQ
        for (k = 0; k < n; k++)
        {
            // Set the number of nonzeros in the kth column of PAQ
            PAQ->p[k] = nz;
            // Column k of PAQ is equal to column S->p[k] of A. j is the
            // starting point for nonzeros and indices for column S->p[k] of A
            j = S->Q_perm[k];
            // Iterate across the nonzeros in column S->p[k]
            for (t = A->p[j]; t < A->p[j+1]; t++)
            {
                // Row i of this nonzero is equal to pinv[A->i[t]]
                PAQ->i[nz] = A->i[t];
                // Move to the next nonzero element of PAQ
                nz++;
            }
        }
    }

    // Finalize the last column of PAQ
    PAQ->p[n] = nz;
    // Set output, return success
    (*PAQ_handle) = PAQ;
    return SPEX_OK;
}
