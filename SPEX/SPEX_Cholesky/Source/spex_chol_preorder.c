//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_preorder: Symbolic ordering and analysis for SPEX Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function performs the symbolic ordering for SPEX Cholesky.
 * Currently, there are three options: user-defined order, COLAMD, or AMD.
 * It is *highly* recommended that AMD is used for Cholesky factorization.
 *
 * Input/output arguments:
 *
 * S:       Symbolic analysis struct. Undefined on input; contains column
 *          permutation and estimates of nnz on output (is the exact
 *          number of nonzeros if AMD is used)
 *
 * A:       Input matrix, unmodified on input/output
 *
 * option:  option->order tells the function which ordering scheme to use
 */

# define SPEX_FREE_ALL                       \
{                                            \
    SPEX_symbolic_analysis_free(&S, option); \
}

#include "spex_chol_internal.h"

SPEX_info spex_chol_preorder
(
    // Output
    SPEX_symbolic_analysis** S_handle,  // Symbolic analysis data structure 
                                        // On input: undefined
                                        // On output: contains the 
                                        // row/column permutation and its
                                        // inverse.
    // Input
    const SPEX_matrix* A,               // Input matrix
    const SPEX_options* option          // Control parameters (use default if NULL)
)
{

    SPEX_info info;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    if ( !spex_initialized() ) return SPEX_PANIC;

    // TODO ensure the following statement is true
    // All inputs have been checked by the caller so ASSERTS are used instead of ifs
    // A can have any data type, but must be in sparse CSC format
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);

    // m = n for Cholesky factorization
    ASSERT(A->n == A->m);

    // Dimension can't be negative
    ASSERT(A->n >= 0); 
    
    // If *S_handle != NULL then it may cause a memory leak.
    ASSERT(*S_handle == NULL);

    //--------------------------------------------------------------------------
    // Allocate symbolic analysis object
    //--------------------------------------------------------------------------
    SPEX_symbolic_analysis* S = NULL;
    // declare indices and dimension of matrix
    int64_t i, k, index, n = A->n;
    int pr = SPEX_OPTION_PRINT_LEVEL(option);

    int64_t anz; // Number of nonzeros in A
    SPEX_matrix_nnz(&anz, A, option);   

    // Allocate memory for S    
    S = (SPEX_symbolic_analysis*) SPEX_calloc(1, sizeof(SPEX_symbolic_analysis));
    if (S == NULL) return SPEX_OUT_OF_MEMORY;

    S->kind = SPEX_CHOLESKY_FACTORIZATION ;

    // Allocate memory for row/column permutation
    S->P_perm = (int64_t*)SPEX_malloc( (n+1)*sizeof(int64_t) );
    if (S->P_perm == NULL)
    {
        SPEX_FREE_ALL;  
        return SPEX_OUT_OF_MEMORY;
    }

    //Check which ordering to use.
    SPEX_col_order order = SPEX_OPTION_ORDER(option);
    switch(order)
    {
        default:
        case SPEX_ORDER_DEFAULT:
        case SPEX_AMD:
        // ---AMD ordering is used (DEFAULT)---
        // S->p is set to AMD's column ordering on A.
        // The number of nonzeros in L is given as AMD's computed
        // number of nonzeros in the Cholesky factor L of A which is the exact
        // nnz(L) for Cholesky factorization (barring numeric cancellation)
        {
            double Control[AMD_CONTROL];           // Declare AMD control
            amd_l_defaults(Control);              // Set AMD defaults
            double Info [AMD_INFO];
            // Perform AMD
            amd_l_order(n, (SuiteSparse_long *)A->p, (SuiteSparse_long *)A->i,
                        (SuiteSparse_long *)S->P_perm, Control, Info);
             S->lnz = Info[AMD_LNZ];        // Exact number of nonzeros for Cholesky
             if (pr > 0)   // Output AMD info if desired
             {
                 SPEX_PRINTF("\n****Column Ordering Information****\n");
                 amd_l_control(Control);
                 amd_l_info(Info);
             }
             double flops=A->n + Info[AMD_NDIV] +2*Info [AMD_NMULTSUBS_LDL]; //n + ndiv + 2*nmultsubs_ldl //CLUSTER
             //printf("%f, %d, ",flops, S->lnz); //CLUSTER
        }
        break;

        case SPEX_NO_ORDERING:
        // ---No ordering is used--- 
        // S->p is set to [0 ... n] and the number of nonzeros in L is estimated 
        // to be 10 times the number of nonzeros in A. 
        // This is a very crude estimate on the nnz(L)
        {
            for (i = 0; i < n+1; i++)
            {
                S->P_perm[i] = i;
            }
            // Very crude estimate for number of L and U nonzeros
            S->lnz = 10*anz;
        }
        break;

        case SPEX_COLAMD:
        // --- COLAMD ordering is used
        // S->p is set as COLAMD's column ordering.
        // The number of nonzeros in L is set as 10 times the number of
        // nonzeros in A. This is a crude estimate.
        {
            // Declared as per COLAMD documentation
            int64_t Alen = 2*anz + 6 *(n+1) + 6*(n+1) + n;
            int64_t* A2 = (int64_t*)SPEX_malloc(Alen*sizeof(int64_t));
            if (!A2)
            {
                // out of memory
                SPEX_FREE_ALL;  
                return SPEX_OUT_OF_MEMORY;
            }
            // Initialize S->p as per COLAMD documentation
            for (i = 0; i < n+1; i++)
            {
                S->P_perm[i] = A->p[i];
            }
            // Initialize A2 per COLAMD documentation
            for (i = 0; i < anz; i++)
            {
                A2[i] = A->i[i];
            }
            int64_t stats[COLAMD_STATS];
            colamd_l(n, n, Alen, (SuiteSparse_long *)A2,
                     (SuiteSparse_long *) S->P_perm, (double *)NULL,
                     (SuiteSparse_long *) stats);
            // estimate for lnz and unz
            S->lnz = 10*anz;

            // Print stats if desired
            if (pr > 0)
            {
                SPEX_PRINTF ("\n****Column Ordering Information****\n");
                colamd_l_report ((SuiteSparse_long *) stats);
                SPEX_PRINTF ("\nEstimated L and U nonzeros: %" PRId64 "\n", S->lnz);
            }
            //Note that A2 is a local-to-this-case variable; so it cannot and should 
            //not be part of the  SPEX_FREE_WORKSPACE or SPEX_FREE_ALL mechanisms
            SPEX_FREE(A2);
        }
        break;

    }

    //--------------------------------------------------------------------------
    // Make sure appropriate space is allocated. It is possible to return
    // estimates which exceed the dimension of L or estimates which are
    // too small for L. In this case, this block of code ensures that the
    // estimates on nnz(L) and nnz(U) are at least n and no more than n*n.
    //--------------------------------------------------------------------------
    // estimate exceeds max number of nnz in A
    if (S->lnz > (double) n*n)
    {
        int64_t nnz = ceil(0.5*n*n);
        S->lnz =  nnz;
    }
    // If estimate < n, it is possible that the first iteration of triangular solve
    // may fail, so we make sure that the estimate is at least n
    if (S->lnz < n)
    {
        S->lnz += n;
    }

    // Allocate pinv
    S->Pinv_perm = (int64_t*)SPEX_calloc(n, sizeof(int64_t));
    if(!(S->Pinv_perm)) return SPEX_OUT_OF_MEMORY;
    
    // Populate pinv
    for (k = 0; k < n; k++)
    {
        index = S->P_perm[k];
        S->Pinv_perm[index] = k;
    }

    //--------------------------------------------------------------------------
    // Set result, report success
    //--------------------------------------------------------------------------
    (*S_handle) = S;
    return SPEX_OK;
}
