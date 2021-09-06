//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Factor: Integer preserving Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

# define SPEX_FREE_ALLOCATION      \
    SPEX_matrix_free(&L, NULL);    \
    SPEX_matrix_free(&rhos, NULL); \

#include "spex_chol_internal.h"  

/* Purpose: This function performs the integer preserving Cholesky factorization.
 * It allows either the left-looking or up-looking integer-preserving Cholesky
 * factorization. In order to compute the L matrix, it performs n iterations of
 * a sparse REF symmetric triangular solve function. The overall factorization
 * is PAP' = LDL'
 * 
 * Importantly, this function assumes that A has already been permuted.
 * 
 * Input arguments of the function:
 * 
 * L_handle:    A handle to the L matrix. Null on input.
 *              On output, contains a pointer to the L matrix.
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. 
 *              On input it contains information that is not used in this 
 *              function such as the row/column permutation
 *              On output it contains the elimination tree and 
 *              the number of nonzeros in L.
 * 
 * rhos_handle: A handle to the sequence of pivots. NULL on input. 
 *              On output it contains a pointer to the pivots matrix.
 *
 * A:           The user's permuted input matrix
 * 
 * left:        A boolean parameter which tells the function whether it is 
 *              performing a left-looking or up-looking factorization. If this
 *              bool is true, a left-looking factorization
 *              is done, otherwise the up-looking factorization is done.
 * 
 * option:      Command options. Notably, option->chol_type indicates whether
 *              it is performing a left-looking (SPEX_CHOL_LEFT) or up-looking 
 *              factorization (SPEX_CHOL_UP) (default)
 * 
 */

//TODO cooments (like up and left) DONE

SPEX_info SPEX_Chol_Factor      
(
    // Output
    SPEX_matrix** L_handle,    // Lower triangular matrix. NULL on input.
    SPEX_matrix** rhos_handle, // Sequence of pivots. NULL on input.
    SPEX_Chol_analysis* S,     // Symbolic analysis struct containing the
                               // elimination tree of A, the column pointers of
                               // L, and the exact number of nonzeros of L.
    // Input
    const SPEX_matrix* A,      // Matrix to be factored   
    const SPEX_options* option //command options
                               // Notably, option->chol_type indicates whether
                               // CHOL_UP (default) or CHOL_LEFT is used.
)
{
    SPEX_info info;
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);
    if (!L_handle || !S || !rhos_handle || !option )
    {
        return SPEX_INCORRECT_INPUT;
    }

    int64_t anz;
    SPEX_CHECK(SPEX_matrix_nnz (&anz, A, option)) ;
    
    if (anz < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare outputs 
    //--------------------------------------------------------------------------

    SPEX_matrix *L = NULL ;
    SPEX_matrix *rhos = NULL ;

    //--------------------------------------------------------------------------
    // Call factorization
    //--------------------------------------------------------------------------

    switch(option->algo) 
    {
        case SPEX_ALGORITHM_DEFAULT:
        case SPEX_CHOL_UP:
            SPEX_CHECK(spex_Chol_Up_Factor(&L, &rhos, S, A, option));
            break;
        case SPEX_CHOL_LEFT:
            SPEX_CHECK(spex_Chol_Left_Factor(&L, &rhos, S, A, option));
            break;
        default:
            return SPEX_INCORRECT_ALGORITHM; //TODO ADD DONE
    }
    /**/
    
    //--------------------------------------------------------------------------
    // Set outputs, return ok
    //--------------------------------------------------------------------------
    (*L_handle) = L;
    (*rhos_handle) = rhos;
    return SPEX_OK;
}
