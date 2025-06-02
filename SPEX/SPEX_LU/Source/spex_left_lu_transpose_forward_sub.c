//------------------------------------------------------------------------------
// SPEX_LU/spex_left_lu_transpose_forward_sub:
//              sparse transpose forward substitution (x = (U'D)\x)
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2025, Christopher Lourenco, Jinhao Chen,,
// Erick Moreno-Centeno, and Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs a sparse transpose roundoff-error-free (REF)
 * forward substitution, that is x = (U'D) \x. This is a subroutine in solving
 * the transposed linear system A'x = b. Mathematically, we do not transpose
 * U directly, instead, since U is stored in CSC format, we can think of it
 * as being U' stored in compressed row format.
 * We also assume that x is dense, thus we do not compute the nonzero pattern
 * and each nonzero in x is iterated across. The system that is solved is
 * thus U' D x_output = x_input, overwriting the right hand side with the
 * solution.
 *
 * On output, the SPEX matrix x is modified.
 */

#define SPEX_FREE_ALL           \
    SPEX_matrix_free(&h, NULL);

#include "spex_lu_internal.h"

SPEX_info spex_left_lu_transpose_forward_sub
(
    const SPEX_matrix U,    // upper triangular matrix
    SPEX_matrix x,          // right hand side matrix of size n*numRHS
    const SPEX_matrix rhos  // sequence of pivots used in factorization
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_REQUIRE(U, SPEX_CSC, SPEX_MPZ);
    SPEX_REQUIRE(x, SPEX_DENSE, SPEX_MPZ);
    SPEX_REQUIRE(rhos, SPEX_DENSE, SPEX_MPZ);

    //--------------------------------------------------------------------------

    int64_t i, hx, k, j, jnew;
    int sgn ;

    // Build the history matrix
    SPEX_matrix h = NULL ;
    SPEX_CHECK (SPEX_matrix_allocate(&h, SPEX_DENSE, SPEX_INT64, x->m, x->n,
        x->nzmax, false, true, NULL));

    // initialize entries of history matrix to be -1
    for (i = 0; i < x->nzmax; i++)
    {
        h->x.int64[i] = -1;
    }


    //--------------------------------------------------------------------------
    // Iterate across each RHS vector
    //--------------------------------------------------------------------------

    for (k = 0; k < x->n; k++)
    {

        //----------------------------------------------------------------------
        // Iterate accross all nonzeros in x. Assume x is dense
        //----------------------------------------------------------------------

        for (i = 0; i < x->m; i++)
        {

            //------------------------------------------------------------------
            // IPGE updates
            //------------------------------------------------------------------

            // Access row i of U'
            // We are finalizing the value of x[i] which is initially set as b[i]
            // Thus, for an example b[i], we are calculating (in the dense case)
            // U[i,1] x[1] + U[i,2] x[2] + ... + U[i, i-1] x[i-1] + U[i,i] x[i] = b[i]
            // Thus U[i,i] x[i] = b[i] - (U[i,1] x[1] + U[i,2] x[2] + ... + U[i, i-1] x[i-1])
            // The first loop iterates through each nonzero in row i of U and performs
            // this submul IPGE update on b[i].
            // Once we have done so, we perform a history update on b[i] to finalize it
            for (j = U->p[i]; j < U->p[i+1]-1; j++)
            {
                // Column index of U[j]
                jnew = U->i[j];
                ASSERT (jnew <= i);

                // Now we history update x[i] if necessary with respect to jnew
                hx = SPEX_2D(h, i, k, int64);
                if (hx < jnew-1)
                {
                    // x[i] = x[i]*rhos[jnew-1]
                    SPEX_MPZ_MUL(SPEX_2D(x, i, k, mpz),
                                 SPEX_2D(x, i, k, mpz),
                                 SPEX_1D(rhos, jnew-1, mpz));

                    if (hx > -1)
                    {
                        SPEX_MPZ_DIVEXACT(SPEX_2D(x,i,k,mpz),
                                          SPEX_2D(x,i,k,mpz),
                                          SPEX_1D(rhos, hx, mpz));
                    }
                }

                // x[i]*rhos[jnew]
                SPEX_MPZ_MUL(SPEX_2D(x,i,k,mpz), SPEX_2D(x,i,k,mpz),
                             SPEX_1D(rhos, jnew, mpz));

                // Now, we update x[i] using U'[i,j]*x[j]
                SPEX_MPZ_SUBMUL(SPEX_2D(x,i,k,mpz),
                                U->x.mpz[j], SPEX_2D(x,jnew,k,mpz));

                if (jnew > 0)
                {
                    // Divide by rhos[jnew-1]
                    SPEX_MPZ_DIVEXACT(SPEX_2D(x,i,k,mpz),
                                  SPEX_2D(x,i,k,mpz),
                                  SPEX_1D(rhos, jnew-1, mpz));
                }
                // Update history
                SPEX_2D(h,i,k,int64) = jnew;
            }
            hx = SPEX_2D(h, i, k, int64);
            // History update to finalize x[i]
            if (hx < i-1)
            {
                // x[i] = x[i] * rhos[i-1]
                SPEX_MPZ_MUL(SPEX_2D(x, i, k, mpz),
                             SPEX_2D(x, i, k, mpz),
                             SPEX_1D(rhos, i-1, mpz));
                // x[i] = x[i] / rhos[hx]
                if (hx > -1)
                {
                    SPEX_MPZ_DIVEXACT(SPEX_2D(x, i, k, mpz),
                                      SPEX_2D(x, i, k, mpz),
                                      SPEX_1D(rhos, hx, mpz));
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Free h memory
    //--------------------------------------------------------------------------

    SPEX_FREE_ALL;
    return SPEX_OK;
}

