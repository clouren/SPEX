//------------------------------------------------------------------------------
// SPEX_Util/SPEX_determine_symmetry: determine if a matrix is symmetric
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Determine if the input A is indeed symmetric prior to
 * factorization.  There are two options as to how to determine the symmetry.
 * By setting the input check_if_numerically_symmetric is true, both the
 * nonzero pattern and the values of the nonzero entries are checked for
 * symmetry. If A passes both of these tests, then we can be sure it is indeed
 * fully symmetric.
 * 
 * If check_if_numerically_symmetric is set to false, only the nonzero pattern
 * of A is checked, thus we cannot gauranteee that the matrix is indeed fully
 * symmetric as the values of the entries is not checked.
 * 
 * If the matrix is determined to be symmetric, SPEX_OK is returned; otherwise,
 * SPEX_UNSYMMETRIC is returned.
 */


// TODO: test coverage

#define SPEX_FREE_ALL               \
{                                   \
    SPEX_matrix_free(&T,NULL);      \
    SPEX_matrix_free(&R,NULL);      \
}

#include "spex_util_internal.h"

SPEX_info SPEX_determine_symmetry
(
    SPEX_matrix* A,
    bool check_if_numerically_symmetric
            // if true, check A=A' (pattern & values). if false,
            // only check if the pattern of A is symmetric, not
            // the values
)
{
    int64_t j;
    SPEX_info info;
    // Declare matrix T
    SPEX_matrix *T = NULL, *R = NULL ;
    // T = A'
    SPEX_CHECK( SPEX_transpose(&T, A));

    // Check if column pointers are the same
    for (j = 0; j <= A->n; j++)
    {
        if (T->p[j] != A->p[j])
        {
            // nnz( A(:,k)) != nnz( A'(:,k))
            SPEX_FREE_ALL ;
            return SPEX_UNSYMMETRIC;
        }
    }

    // R = T'
    SPEX_CHECK( SPEX_transpose(&R, T));
    // then compare R and T

    // Check if i values are the same
    for (j = 0; j < R->nz; j++)
    {
        if (T->i[j] != R->i[j])
        {
            // A[i][j] != A[j][i], unsymmetric
            SPEX_FREE_ALL ;
            return SPEX_UNSYMMETRIC;
        }
    }

    // If we are performing an exhaustive search, we check the x values as well
    // This is by far the most expensive part of checking the symmetry.
    if (check_if_numerically_symmetric)
    {
        int r;
        for (j = 0; j < R->nz; j++)
        {
            SPEX_CHECK(SPEX_mpz_cmp(&r, R->x.mpz[j], T->x.mpz[j]));
            if ( r != 0)
            {
                // Pattern is symmetric, values are not
                SPEX_FREE_ALL ;
                return SPEX_UNSYMMETRIC;
            }
        }
    }
    SPEX_FREE_ALL ;
    return SPEX_OK;
}

