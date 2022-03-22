//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Factor: Integer preserving Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

# define SPEX_FREE_ALL                   \
{                                        \
    SPEX_factorization_free(&F, option); \
}


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
 * F_handle:    A handle to the factorization struct. Null on input.
 *              On output, contains a pointer to the factorization (this 
 *              includes matrix L)
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. 
 *              On input it contains the elimination tree and 
 *              the number of nonzeros in L.
 *
 * A:           The user's permuted input matrix
 * 
 * option:      Command options. Notably, option->chol_type indicates whether
 *              it is performing a left-looking (SPEX_CHOL_LEFT) or up-looking 
 *              factorization (SPEX_CHOL_UP) (default)
 * 
 */


SPEX_info spex_chol_factor      
(
    // Output
    SPEX_factorization **F_handle, // Cholesky factorization
    //Input
    const SPEX_symbolic_analysis* S, // Symbolic analysis struct containing the
                               // elimination tree of A, the column pointers of
                               // L, and the exact number of nonzeros of L.
    const SPEX_matrix* A,      // Matrix to be factored   
    const SPEX_options* option //command options
                               // Notably, option->chol_type indicates whether
                               // CHOL_UP (default) or CHOL_LEFT is used.
)
{
    SPEX_info info;

    SPEX_factorization *F = NULL ;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);
    /**/if (!F_handle || !S || !option )
    {
        return SPEX_INCORRECT_INPUT;
    }

    int64_t anz;
    SPEX_CHECK(SPEX_matrix_nnz (&anz, A, option)) ;
    
    if (anz < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }/**/ //TOCHECK why do I need to comment this so that things will compile

    (*F_handle) = NULL ;

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    int64_t n = A->n ;

    //SPEX_factorization *F = NULL ;
    // allocate memory space for the factorization
    F = (SPEX_factorization*) SPEX_calloc(1, sizeof(SPEX_factorization));
    if (F == NULL)
    {
        return SPEX_OUT_OF_MEMORY;
    }
    SPEX_CHECK(SPEX_factorization_create(&F,option));
    // set factorization kind
    F->kind = SPEX_CHOLESKY_FACTORIZATION;

    // Allocate and set scale_for_A
    SPEX_CHECK(SPEX_mpq_init(F->scale_for_A));
    SPEX_CHECK(SPEX_mpq_set (F->scale_for_A, A->scale));

    // Inverse pivot ordering
    F->Pinv_perm = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
    // Actual row permutation, the inverse of pinv. This
    // is used for sorting
    F->P_perm =    (int64_t*) SPEX_malloc (n * sizeof(int64_t));
    // column permutation, to be copied from S->Q_perm
    F->Q_perm =    (int64_t*) SPEX_malloc (n * sizeof(int64_t));


    // copy column permutation from symbolic analysis to factorization
    memcpy(F->Q_perm, S->Q_perm, n * sizeof(int64_t));


    //--------------------------------------------------------------------------
    // Declare memory for rhos, L, and U
    //--------------------------------------------------------------------------
    // Create rhos, a global dense mpz_t matrix of dimension n*1
    SPEX_CHECK (SPEX_matrix_allocate(&(F->rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, false, option));

    // Allocate L without initializing each entry.
    // L is allocated to have nnz(L) which is estimated by the symbolic
    // analysis. However, unlike traditional matrix allocation, the second
    // boolean parameter here is set to false, so the individual values of
    // L are not allocated. Instead, a more efficient method to
    // allocate these values is done in the factorization to reduce
    // memory usage.
    SPEX_CHECK (SPEX_matrix_allocate(&(F->L), SPEX_CSC, SPEX_MPZ, n, n, S->lnz,
        false, false, option));

    //--------------------------------------------------------------------------
    // Call factorization
    //--------------------------------------------------------------------------

    switch(option->algo) 
    {
        case SPEX_ALGORITHM_DEFAULT:
        case SPEX_CHOL_UP:
            SPEX_CHECK(spex_chol_up_factor(&(F->L), &(F->rhos), S, A, option));
            break;
        case SPEX_CHOL_LEFT:
            SPEX_CHECK(spex_chol_left_factor(&(F->L), &(F->rhos), S, A, option));
            break;
        default:
            return SPEX_INCORRECT_ALGORITHM; 
    }
    /**/
    memcpy(F->Pinv_perm, S->Pinv_perm, n * sizeof(int64_t));

    //--------------------------------------------------------------------------
    // Set outputs, return ok
    //--------------------------------------------------------------------------
    (*F_handle) = F ;
    return SPEX_OK;
}
