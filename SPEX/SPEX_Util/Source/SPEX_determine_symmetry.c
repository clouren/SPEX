///------------------------------------------------------------------------------
// SPEX_Util/SPEX_determine_symmetry: This function determines if the input matrix is symmetric.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Determine if the input A is indeed symmetric prior to factorization.
 * There are two options as to how to determine the symmetry. 
 * By setting the input exhaustive true, both the nonzero pattern and the values
 * of the nonzero entries are checked for symmetry. If A passes both of these tests,
 * then we can be sure it is indeed fully symmetric.
 * 
 * If exhaustive is set to false, only the nonzero pattern of A is checked,
 * thus we cannot gauranteee that the matrix is indeed fully symmetric as the values
 * of the entries is not checked.
 * 
 * If the matrix is determined to be symmetric, true is returned; otherwise, false is
 * returned.
 * 
 */

#include "spex_util_internal.h"
//TODO fix me
bool SPEX_determine_symmetry
(
    SPEX_matrix* A,
    bool exhaustive
)
{
    int64_t j;
    
    // Declare matrix T
    SPEX_matrix *T = NULL;    
    // T = A'
    SPEX_transpose(&T, A);
    
    // Check if i values are the same
    for (j = 0; j < A->nz; j++)
    {
        if (T->i[j] != A->i[j])
        {
            //printf("\nError, matrix is not symmetric\n");
            SPEX_matrix_free(&T,NULL);
            return false;
        }
    }
    
    // Check if column pointers are the same
    for (j = 0; j <= A->n; j++)
    {
        if (T->p[j] != A->p[j])
        {
            //printf("\nError, matrix is not symmetric\n");
            SPEX_matrix_free(&T,NULL);
            return false;
        }
    }

    // If we are performing an exhaustive search, we check the x values as well
    // This is by far the most expensive part of checking the symmetry.
    if (exhaustive == true)
    {
        int r;
        for (j = 0; j < A->nz; j++)
        {
            SPEX_mpz_cmp(&r, A->x.mpz[j], T->x.mpz[j]);
            if ( r != 0)
            {
                //printf("\nError, pattern is symmetric, values are not\n");
                SPEX_matrix_free(&T,NULL);
                return false;
            }
        }
    }
    SPEX_matrix_free(&T,NULL);
    return true;
        
}
