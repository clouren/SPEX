//------------------------------------------------------------------------------
// SPEX_QR/SPEX_qr_analyze: Perform the symbolic analysis for QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Lorena Mejia Domenzain, Christopher Lourenco,
// Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: perform the symbolic analysis of A for the QR factorization,
 * that is, preordering A, computing the elimination tree, getting the column
 * counts of ATA, setting the column pointers and number of non zeros of R.
 *
 * Input arguments of the function:
 *
 * S:           Symbolic analysis struct for QR factorization.
 *              On input it's NULL
 *              On output it contains the row/column permutation, the elimination
 *              tree, and the number of nonzeros in R.
 *
 * A:           User's input matrix (Must be SPEX_MPZ and SPEX_CSC)
 *
 * option:      Command options (Default if NULL)
 *
 */

#define SPEX_FREE_WORKSPACE         \
{                                   \
    SPEX_matrix_free(&PAQ, NULL);   \
    SPEX_FREE(c);                   \
    SPEX_FREE(cInv);                \
}

#define SPEX_FREE_ALL                               \
{                                                   \
    SPEX_FREE_WORKSPACE ;                           \
    SPEX_symbolic_analysis_free (&S, option);      \
}

#include "spex_cholesky_internal.h"
#include "spex_qr_internal.h"

SPEX_info SPEX_qr_analyze
(
    // Output
    SPEX_symbolic_analysis *S_handle, // Symbolic analysis data structure
    // Input
    const SPEX_matrix A,        // Input matrix. Must be SPEX_MPZ and SPEX_CSC
    const SPEX_options option   // Command options (Default if NULL)
)
{

    SPEX_info info;
    // SPEX must be initialized
    if (!spex_initialized())
    {
        return SPEX_PANIC;
    }

    // Check inputs
    if ( !S_handle || !A)
    {
        printf("hcin\n");
        return SPEX_INCORRECT_INPUT;
    }

    // SPEX must be CSC
    SPEX_REQUIRE_KIND(A, SPEX_CSC);

    // Declare permuted matrix and S
    SPEX_matrix PAQ = NULL;
    SPEX_symbolic_analysis S = NULL;
    // Declare local variables for symbolic analysis
    int64_t n = A->n, m=A->m;
    int64_t *post = NULL;
    int64_t *c = NULL, *cInv=NULL;
    int64_t i, nz;

    //--------------------------------------------------------------------------
    // Preorder: obtain the row/column ordering of ATA (Default is COLAMD)
    //--------------------------------------------------------------------------

    SPEX_CHECK( spex_qr_preorder(&S, A, option) );

    //--------------------------------------------------------------------------
    // Permute matrix A, that is apply the row/column ordering from the
    // symbolic analysis step to get the permuted matrix PAQ.
    //--------------------------------------------------------------------------

    SPEX_CHECK( spex_qr_permute_A(&PAQ, A, true, S, option) ); //TODO can make false when you can transpose an empty matrix 
    //SPEX_matrix_check(PAQ, option);
    
    //--------------------------------------------------------------------------
    // Symbolic Analysis: compute the elimination tree of PAQ
    //--------------------------------------------------------------------------

    // Obtain elimination tree of A
    SPEX_CHECK( spex_qr_etree(&S->parent, PAQ) );
    

    // Postorder the elimination tree of A
    SPEX_CHECK( spex_cholesky_post(&post, S->parent, n) );


    // Get the column counts of A
    SPEX_CHECK( spex_qr_counts(&c, PAQ, S->parent, post) ); //c is S->cp but backwards

    // Set the column pointers of R
    S->cp = (int64_t*) SPEX_malloc( (n+1)*sizeof(int64_t*));
    if (S->cp == NULL)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    cInv = (int64_t*) SPEX_malloc(n* sizeof (int64_t));
    for(i=1;i<=n;i++)
    {
        cInv[i]=c[n-i];//FIXME there has to be a better way of doing this check L vs R
    }
    
    for(i=0;i<=n;i++)
    {
        //cInv[i]=c[n-i];//FIXME there has to be a better way of doing this check L vs R
        printf("%ld cinv %ld c %ld\n",i, cInv[i],c[i]);
    }
    
    //FIXME hardcode col counts for miniZeros, something is wrong with the analysis
    /*cInv[1]=1;
    cInv[2]=1;
    cInv[3]=2;
    cInv[4]=2;
    cInv[5]=3;*/

    SPEX_CHECK( spex_cumsum(S->cp, cInv, n));

    nz=0;
    for (i = 0 ; i <= n ; i++)
    {
        nz += cInv[i];//S->cp [i] ;
    }
    S->unz=nz;//suma de todos los elementos de c

    //TODO non zero patern of Q?? dense??

    // set num non-zeros in Q
    S->lnz = m*n; //Q is dense right now

    //--------------------------------------------------------------------------
    // Set output, free all workspace and return success
    //--------------------------------------------------------------------------

    (*S_handle) = S ;
    SPEX_FREE_WORKSPACE ;
    return (SPEX_OK);
}
