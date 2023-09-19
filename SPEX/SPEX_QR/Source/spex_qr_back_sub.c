//------------------------------------------------------------------------------
// SPEX_QR/Source/spex_qr_basic_solve.c: Basic solution back solve
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF backward substitution for 
 * underdetermined SLEs, solving the system Rx = b. Where the last n-rank rows of
 * R are 0.
 *
 * R is a sparse mpz matrix, and bx is a dense mpz matrix.  The diagonal entry
 * of U must appear as the last entry in each column.
 *
 * The input argument bx contains b on input, and it is overwritten on output
 * by the solution x.
 */

# include "spex_qr_internal.h"

SPEX_info spex_qr_back_sub  // performs sparse REF backward substitution
(
    SPEX_matrix bx,         // right hand side matrix
    const SPEX_matrix R,   // input upper triangular matrix
    const int64_t rank     // rank of right triangular matrix
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_REQUIRE (R,  SPEX_CSC,   SPEX_MPZ);
    SPEX_REQUIRE (bx, SPEX_DENSE, SPEX_MPZ);

    //--------------------------------------------------------------------------

    int sgn;
    mpz_t *Rx = R->x.mpz;
    int64_t *Ri = R->i;
    int64_t *Rp = R->p;

    for (int64_t k = 0; k < bx->n; k++)
    {
        // Start at bx[n]
        for (int64_t j = rank-1; j >= 0; j--)
        {
            // If bx[j] is zero skip this iteration
            SPEX_MPZ_SGN(&sgn, SPEX_2D( bx, j, k, mpz));
            if (sgn == 0) {continue;}

            // Obtain bx[j]
            SPEX_MPZ_DIVEXACT(SPEX_2D(bx, j, k, mpz),
                                          SPEX_2D(bx, j, k, mpz),
                                          Rx[Rp[j+1]-1]);
            for (int64_t i = Rp[j]; i < Rp[j+1]-1; i++)
            {
                SPEX_MPZ_SGN(&sgn, Rx[i]);
                if (sgn == 0) {continue;}
                // bx[i] = bx[i] - Rx[i]*bx[j]
                SPEX_MPZ_SUBMUL(SPEX_2D(bx, Ri[i], k, mpz),
                                Rx[i], SPEX_2D(bx, j, k, mpz));
            }
        }
    }

    return (SPEX_OK);
}
