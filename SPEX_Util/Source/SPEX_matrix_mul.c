//------------------------------------------------------------------------------
// SPEX_Util/spex_matrix_mul: multiplies a matrix by a scalar
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function multiplies the matrix x (CSC, triplet, or dense) by a
 * scalar.  This function requires x to have type mpz_t.
 *
 * On output the values of x are modified.
 */

#include "spex_util_internal.h"

SPEX_info SPEX_matrix_mul   // multiplies x by a scalar
(
    SPEX_matrix *x,         // matrix to be multiplied
    const mpz_t scalar      // scalar to multiply by
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_REQUIRE_TYPE (x, SPEX_MPZ) ;

    //--------------------------------------------------------------------------
    // x = x * scalar
    //--------------------------------------------------------------------------

    int64_t nz = SPEX_matrix_nnz (x, NULL) ;
    for (int64_t i = 0; i < nz; i++)
    {
        // x[i] = x[i]*scalar
        SPEX_CHECK( SPEX_mpz_mul( x->x.mpz[i], x->x.mpz[i], scalar));
    }

    return (SPEX_OK) ;
}

