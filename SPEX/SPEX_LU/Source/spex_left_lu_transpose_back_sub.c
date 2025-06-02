//------------------------------------------------------------------------------
// SPEX_LU/spex_left_lu_transpose_back_sub:
//              sparse transpose REF backward substitution (x = L'\x)
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2025, Christopher Lourenco, Jinhao Chen,,
// Erick Moreno-Centeno, and Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse transposed REF backward substitution, solving
 * the system L'x = b. This is a subroutine in the transpose solve A' x = b
 *
 * Note that prior to this, x is multiplied by
 * the determinant of A. Thus a standard substitution can be used.
 *
 * L is a sparse mpz matrix, and bx is a dense mpz matrix.
 *
 * The input argument bx contains b on input, and it is overwritten on output
 * by the solution x.
 */


#include "spex_lu_internal.h"

SPEX_info spex_left_lu_transpose_back_sub  // performs sparse transpose REF backward sub
(
    const SPEX_matrix L,    // input lower triangular matrix
    const SPEX_matrix rhos, // Sequence of pivots
    SPEX_matrix bx          // right hand side matrix
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_REQUIRE (L,  SPEX_CSC,   SPEX_MPZ);
    SPEX_REQUIRE (bx, SPEX_DENSE, SPEX_MPZ);

    //--------------------------------------------------------------------------

    int sgn;
    mpz_t *Lx = L->x.mpz;
    int64_t *Li = L->i;
    int64_t *Lp = L->p;

    for (int64_t k = 0; k < bx->n; k++)
    {
        // Start at bx[n]
        for (int64_t j = L->n-1; j >= 0; j--)
        {

            for (int64_t i = Lp[j]; i < Lp[j+1]; i++)
            {
                // Since row indices are not sorted, ensure
                // that we only use entries that have a row
                // index greater than j
                if (Li[i] > j)
                {
                    SPEX_MPZ_SUBMUL(SPEX_2D(bx, j, k, mpz),
                                Lx[i], SPEX_2D(bx, Li[i], k, mpz));
                }
            }
            SPEX_MPZ_DIVEXACT(SPEX_2D(bx,j,k,mpz),
                              SPEX_2D(bx,j,k,mpz),
                              rhos->x.mpz[j]);
        }
    }

    return (SPEX_OK);
}

