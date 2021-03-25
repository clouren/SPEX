//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_analyze: symbolic ordering and analysis for sparse Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs the symbolic ordering for SPEX Cholesky.
 * Currently, there are three options: user-defined order, COLAMD, or AMD.
 * It is *highly* recommended that AMD is used for Cholesky factorization.
 *
 * Input/output arguments:
 *
 * S:       Symbolic analysis struct. Undefined on input; contains column
 *          permutation and estimates of nnz(L) and nnz(U) nnz on output
 *
 * A:       Input matrix, unmodified on input/output
 *
 * option:  option->order tells the function which ordering scheme to use
 *
 */

// SPEX_Chol_analyze creates the SPEX_Chol_analysis object S_Chol.  Use
// SPEX_Chol_analysis_free to delete it.
// This function is a slightly modified version of SPEX_Left_LU's 
// LU analysis function

#include "spex_chol_internal.h"

SPEX_info SPEX_Chol_preorder
(
    SPEX_Chol_analysis** S_handle, // symbolic analysis (row/column perm. and nnz L,U)
    const SPEX_matrix *A,        // Input matrix
    const SPEX_options *option   // Control parameters, if NULL, use default
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    // A can have any data type, but must be in sparse CSC format
    SPEX_REQUIRE_KIND (A, SPEX_CSC) ;

    // m = n for Cholesky factorization
    ASSERT( A->n == A->m);
    
    if (!S_handle || A->n != A->m)
    {
        return SPEX_INCORRECT_INPUT;
    }
    (*S_handle) = NULL ;

    //--------------------------------------------------------------------------
    // allocate symbolic analysis object
    //--------------------------------------------------------------------------

    SPEX_Chol_analysis *S = NULL ;
    int64_t i, n = A->n;
    ASSERT( n >= 0); // Dimension can't be negative
    int64_t anz;
    SPEX_matrix_nnz(&anz, A, option);   // Number of nonzeros in A
    
    // Allocate memory for S
    S = (SPEX_Chol_analysis*) SPEX_malloc(sizeof(SPEX_Chol_analysis));
    if (S == NULL) {return SPEX_OUT_OF_MEMORY;}

    // Allocate memory for column permutation
    S->q = (int64_t*) SPEX_malloc((n+1) * sizeof(int64_t));
    if (S->q == NULL)
    {
        SPEX_FREE(S);
        return SPEX_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // No ordering is used. S->q is set to [0 ... n] and the number of nonzeros
    // in L is estimated to be 10 times the number of nonzeros in A. This
    // is a very crude estimate on the nnz(L)
    //--------------------------------------------------------------------------

    SPEX_col_order order = SPEX_OPTION_ORDER (option) ;
    int pr = SPEX_OPTION_PRINT_LEVEL (option) ;

    if (order == SPEX_NO_ORDERING)
    {
        for (i = 0; i < n+1; i++)
        {
            S->q[i] = i;
        }
        // estimates for number of L and U nonzeros
        S->lnz = 10*anz;
    }

    //--------------------------------------------------------------------------
    // The COLAMD ordering is used. S->q is set as COLAMD's column ordering.
    // The number of nonzeros in L is set as 10 times the number of
    // nonzeros in A. This is a crude estimate.
    //--------------------------------------------------------------------------
    else if (order == SPEX_AMD)
    {
        // Declared as per COLAMD documentation
        int64_t Alen = 2*anz + 6 *(n+1) + 6*(n+1) + n;
        int64_t* A2 = (int64_t*) SPEX_malloc(Alen* sizeof(int64_t));
        if (!A2)
        {
            // out of memory
            SPEX_Chol_analysis_free (&S) ;
            return (SPEX_OUT_OF_MEMORY) ;
        }
        // Initialize S->q as per COLAMD documentation
        for (i = 0; i < n+1; i++)
        {
            S->q[i] = A->p[i];
        }
        // Initialize A2 per COLAMD documentation
        for (i = 0; i < anz; i++)
        {
            A2[i] = A->i[i];
        }
        int64_t stats [COLAMD_STATS];
        colamd_l (n, n, Alen, (SuiteSparse_long *) A2,
            (SuiteSparse_long *) S->q, (double *) NULL,
            (SuiteSparse_long *) stats) ;
        // estimate for lnz and unz
        S->lnz = 10*anz;

        // Print stats if desired
        if (pr > 0)
        {
            SPEX_PRINTF ("\n****Column Ordering Information****\n") ;
            colamd_l_report ((SuiteSparse_long *) stats) ;
            SPEX_PRINTF ("\nEstimated L and U nonzeros: %" PRId64 "\n", S->lnz);
        }
        SPEX_FREE(A2);
    }

    //--------------------------------------------------------------------------
    // The AMD ordering is used. S->q is set to AMD's column ordering on
    // A. The number of nonzeros in L is given as AMD's computed
    // number of nonzeros in the Cholesky factor L of A which is the exact
    // nnz(L) for Cholesky factorization (barring numeric cancellation)
    //--------------------------------------------------------------------------
    else
    {
        double Control [AMD_CONTROL];           // Declare AMD control
        amd_l_defaults (Control) ;              // Set AMD defaults
        double Info [AMD_INFO];
        // Perform AMD
        amd_l_order(n, (SuiteSparse_long *) A->p, (SuiteSparse_long *) A->i,
            (SuiteSparse_long *) S->q, Control, Info) ;
        S->lnz = Info[AMD_LNZ];        // estimate for unz and lnz
        if (pr > 0)   // Output AMD info if desired
        {
            SPEX_PRINTF ("\n****Column Ordering Information****\n") ;
            amd_l_control (Control) ;
            amd_l_info (Info) ;
        }
    }

    //--------------------------------------------------------------------------
    // Make sure appropriate space is allocated. It's possible to return
    // estimates which exceed the dimension of L and U or estimates which are
    // too small for L U. In this case, this block of code ensures that the
    // estimates on nnz(L) and nnz(U) are at least n and no more than n*n.
    //--------------------------------------------------------------------------
    // estimate exceeds max number of nnz in A
    if (S->lnz > (double) n*n)
    {
        int64_t nnz = ceil(0.5*n*n);
        S->lnz =  nnz;
    }
    // If estimate < n, first iteration of left-looking TS may fail
    if (S->lnz < n)
    {
        S->lnz = S->lnz + n;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*S_handle) = S ;
    return SPEX_OK;
}
