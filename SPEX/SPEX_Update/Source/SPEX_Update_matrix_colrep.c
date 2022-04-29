//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_matrix_colrep: perform a single column replacement
// for target matrix
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function swaps a specified column of a given m-by-n matrix
 * with the column from a m-by-1 matrix. Both matrices must be of
 * SPEX_DYNAMIC_CSC SPEX_MPZ. On ouput, both matrices are modified.
 */


#include "spex_update_internal.h"

SPEX_info SPEX_Update_matrix_colrep// performs column replacement
(
    SPEX_matrix *A,         // m-by-n target matrix of SPEX_DYNAMIC_CSC MPZ
    SPEX_matrix *vk,        // m-by-1 SPEX_DYNAMIC_CSC MPZ matrix that contains
                            // the column vector to replace the k-th column of A
    int64_t k,              // The column index that vk will be inserted, 0<=k<n
    const SPEX_options *option// Command parameters
)
{
    SPEX_info info ;

    SPEX_REQUIRE(A , SPEX_DYNAMIC_CSC, SPEX_MPZ);
    SPEX_REQUIRE(vk, SPEX_DYNAMIC_CSC, SPEX_MPZ);
    if (A->m != vk->m || k < 0 || k >= A->n || vk->n != 1)
    {
        return SPEX_INCORRECT_INPUT;
    }

    SPEX_vector *Vtmp = A->v[k];
    A->v[k] = vk->v[0];
    vk->v[0] = Vtmp;
    return SPEX_OK;
}
