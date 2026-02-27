//------------------------------------------------------------------------------
// SPEX_LU/SPEX_lu_rank: compute the rank of sparse square A
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2026, Christopher Lourenco, Jinhao Chen,,
// Erick Moreno-Centeno, and Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This code utilizes the SPEX Left LU factorization to exactly
 * compute the rank of A. A must be square
 *
 * Input/Output arguments:
 *
 * x_handle:    A pointer to the solution of the linear system. The output is
 *              allowed to be returned in either double precision, mpfr_t, or
 *              rational mpq_t
 *
 * type:        Data structure of output desired. Must be either SPEX_MPQ,
 *              SPEX_FP64, or SPEX_MPFR
 *
 * A:           User's input matrix. It must be populated prior to calling this
 *              function.
 *
 * b:           Collection of right hand side vectors. Must be populated prior
 *              to factorization.
 *
 * option:      Struct containing various command parameters for the
 *              factorization. If NULL on input, default values are used.
 */

# define SPEX_FREE_WORKSPACE                    \
    SPEX_symbolic_analysis_free (&S, option);

# define SPEX_FREE_ALL              \
    SPEX_FREE_WORKSPACE             \

#include "spex_lu_internal.h"

SPEX_info SPEX_lu_rank
(
        // Output
        int64_t* rank,                  // rank of A
        // Input
        const SPEX_matrix A,            // Input matrix
        const SPEX_options option       // Command options
)
{

    //-------------------------------------------------------------------------
    // check inputs
    //-------------------------------------------------------------------------

    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC);
    
    // get option->algo, or use SPEX_ALGORITHM_DEFAULT if option is NULL:
    SPEX_factorization_algorithm algo = SPEX_OPTION_ALGORITHM(option);
    if (algo != SPEX_ALGORITHM_DEFAULT && algo != SPEX_LU_LEFT)
    {
        return SPEX_INCORRECT_ALGORITHM;
    }

    SPEX_REQUIRE (A, SPEX_CSC,   SPEX_MPZ);

    SPEX_symbolic_analysis S = NULL;
    SPEX_factorization F = NULL ;

    //--------------------------------------------------------------------------
    // Symbolic Analysis
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_lu_analyze(&S, A, option));

    //--------------------------------------------------------------------------
    // LU Factorization to find rank
    //--------------------------------------------------------------------------

    SPEX_CHECK(spex_lu_rank_factorize(rank, A, S, option));
    
    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SPEX_FREE_WORKSPACE ;
    return (SPEX_OK);
}

