//------------------------------------------------------------------------------
// SPEX_Util/spex_matrix_div: divide a matrix by a scalar
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function takes as input a dense SPEX_matrix x, which is MPZ,
 * and divides it a scalar.  This division is then stored in a dense MPQ
 * matrix, which must have the same number of entries as x.  This is used
 * internally to divide the solution vector by the determinant of the matrix.
 *
 * On output, the contents of the matrix x2 are modified.
 */

#define SPEX_FREE_WORK            \
    SPEX_MPQ_CLEAR(scalar2);

#define SPEX_FREE_ALL               \
    SPEX_FREE_WORK                  \
    SPEX_matrix_free (&x2, NULL) ;

#include "spex_util_internal.h"

SPEX_info SPEX_matrix_div // divides the x matrix by a scalar
(
    SPEX_matrix **x2_handle,    // x2 = x/scalar
    SPEX_matrix* x,             // input vector x
    const mpz_t scalar,         // the scalar
    const SPEX_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_matrix *x2 = NULL ;
    (*x2_handle) = NULL ;
    SPEX_REQUIRE (x, SPEX_DENSE, SPEX_MPZ) ;

    //--------------------------------------------------------------------------
    // Set scalar2 = scalar
    //--------------------------------------------------------------------------

    mpq_t scalar2 ;
    SPEX_MPQ_SET_NULL (scalar2) ;
    SPEX_CHECK (SPEX_mpq_init (scalar2)) ;
    SPEX_CHECK (SPEX_mpq_set_num (scalar2, scalar)) ;

    //--------------------------------------------------------------------------
    // allocate x2
    //--------------------------------------------------------------------------

    SPEX_CHECK (SPEX_matrix_allocate(&x2, SPEX_DENSE, SPEX_MPQ, x->m, x->n,
        0, false, true, option)) ;

    //--------------------------------------------------------------------------
    // iterate each entry of x, copy to x2 and divide it by scalar
    //--------------------------------------------------------------------------

    int64_t nz = SPEX_matrix_nnz (x, NULL) ;
    for (int64_t i = 0; i < nz; i++)
    {
        // Set x2[i] = x[i]
        SPEX_CHECK (SPEX_mpq_set_num (x2->x.mpq[i], x->x.mpz[i])) ;
        // x2[i] = x2[i] / scalar2
        SPEX_CHECK (SPEX_mpq_div (x2->x.mpq[i], x2->x.mpq[i], scalar2)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*x2_handle) = x2 ;
    return (SPEX_OK) ;
}

