//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR_dense.c: Dense REF QR factorization
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* This code contains a dense REF QR factorization. It is meant to be proof of
 * concept for the REF QR factorization algorithms. It computes the
 * factorization A = Q D R where Q and R is integer. Note that all code in this
 * version of SPEX QR assumes that A is a fully dense matrix; thus these
 * routines are not yet appropriate for sparse matrices.
 */


# include "SPEX.h"
# include "spex_util_internal.h"
#include "spex_demos.h"

#define FREE_WORKSPACE                          \
    SPEX_matrix_free(&A,  NULL);                \
    SPEX_matrix_free(&A2, NULL);                \
    SPEX_matrix_free(&R,  NULL);                \
    SPEX_matrix_free(&Q,  NULL);                \
    SPEX_matrix_free(&R2, NULL);                \
    SPEX_matrix_free(&Q2, NULL);                \
    SPEX_FREE(option);                          \
    SPEX_finalize();                            \

#ifndef ASSERT
#define ASSERT assert
#endif

int main( int argc, char *argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX QR, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_initialize();

    // Input arguments.
    unsigned int seed;
    int64_t m;
    int64_t n;
    int64_t lower;
    int64_t upper;

    // Process the command line. If there are not exactly 5 input arguments,
    // then default values are used.
    if (argc != 6)
    {
        printf("\nExpected usage: ./SPEX_QR_dense SEED M N LOWER UPPER\n");
        printf("\nUsing default settings\n");
        seed = 10;
        m = 100;
        n = 50;
        lower = 1;
        upper = 10;
    }
    // Acquire input arguments
    else
    {
        seed = (unsigned int) atoi(argv[1]);
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        lower = atoi(argv[4]);
        upper = atoi(argv[5]);
    }

    m=4;n=4;seed=14;
    // Input checks
    ASSERT(m >= 0);
    ASSERT(n >= 0);
    ASSERT(lower < upper);
    ASSERT(seed >= 1);

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    SPEX_info ok;
    SPEX_matrix A = NULL ;     // Integer matrix to be factorized
    SPEX_matrix A2 = NULL;     // Matrix to be randomly generated

    SPEX_matrix Ainit = NULL;     // Matrix to be randomly generated

    // Next we define 3 Q R pairs. Each pair is generated via a different dense
    // algorithm
    SPEX_matrix Q = NULL;
    SPEX_matrix R = NULL;
    SPEX_matrix Q2 = NULL;
    SPEX_matrix R2 = NULL;
    
    // SPEX Options
    SPEX_options option = NULL;
    SPEX_create_default_options(&option);
    if (!option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Generate a random dense matrix
    //--------------------------------------------------------------------------

    SPEX_generate_random_matrix ( &Ainit, m, n, seed, lower, upper);
    Ainit->nz = m*n;

    // Create A as a copy of Ainit
    // A is a copy of the Ainit matrix. A is a sparse matrix with mpz_t entries
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, Ainit, option);

    // Create A2 as a copy of Ainit
    // A2 is a copy of the Ainit matrix. A is a dense matrix with mpz_t entries
    SPEX_matrix_copy(&A2, SPEX_DENSE, SPEX_MPZ, Ainit, option);

    //--------------------------------------------------------------------------
    // Factorize
    //--------------------------------------------------------------------------


    SPEX_QR_IPGE( A2, &R2, &Q2);

    // playing with sparse
    //will need to change DEMO_OK back to SPEX_CHECK
    SPEX_info info;
    SPEX_matrix rhos = NULL,R3=NULL;
    int64_t *h;
    int64_t j=0, nnz=n*m;
    int64_t i,pQ, pR;
    int sgn;

    h = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    DEMO_OK (SPEX_matrix_allocate(&(rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));

    // Allocate R. We are performing the Thin REF QR factorization so
    // R is n*n
    //DEMO_OK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, n*n, false, true, NULL));
    SPEX_matrix_copy(&R, SPEX_CSC, SPEX_MPZ, R2, option);

    //Q=A (will need to change as soon as A is actually sparse)
    DEMO_OK(SPEX_matrix_copy(&Q, SPEX_CSC, SPEX_MPZ, A, NULL));
    //will also need to change once A is actually sparse (and we're actually only doing one col)
    for (i = 0; i < nnz; i++)
    {
        h[i] = -1;
    }

    for (i = 0; i < n; i++)
    {
        SPEX_MPZ_SET(rhos->x.mpz[i], SPEX_2D(R2, i, i, mpz)); //rho^i=R2(i,i)
    }

    option->print_level = 3;
    printf("printA: \n");
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, Ainit, option);
    SPEX_matrix_check(A, option); //checa que onda, a lo mejor es lo de option
    printf("finishA\n");

    
    //R

   

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;

}

