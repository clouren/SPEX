//------------------------------------------------------------------------------
// Demos/SPEX_Update_demo_chol.c: demo for SPEX_Update library
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* A simple example to show how to perform rank-1 Cholesky factorization
 * update/downdate.
 */

#include "demos.h"

#define FREE_WORKSPACE                           \
{                                                \
    SPEX_matrix_free(&A, option);                \
    SPEX_matrix_free(&w, option);                \
    SPEX_symbolic_analysis_free(&S, option);     \
    SPEX_factorization_free(&F, option);         \
    SPEX_FREE(option);                           \
    SPEX_finalize();                             \
}

int main( int argc, char *argv[] )
{
    SPEX_info ok;
    //--------------------------------------------------------------------------
    // Initialize SPEX CHOLMOD process
    //--------------------------------------------------------------------------

    SPEX_initialize () ;

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    SPEX_options option = NULL;
    SPEX_factorization F = NULL;
    SPEX_symbolic_analysis S = NULL;
    SPEX_matrix A = NULL;
    SPEX_matrix w = NULL;
    DEMO_OK(SPEX_create_default_options(&option));

    //--------------------------------------------------------------------------
    // read matrix and store as a SPEX_CSC SPEX_MPZ matrix A
    //--------------------------------------------------------------------------
    char *mat_name, *rhs_name;
    int64_t rat = 1;
    DEMO_OK(SPEX_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));
    
    FILE *mat_file = fopen(mat_name, "r");
    if (mat_file == NULL)
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    DEMO_OK(SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);

    //--------------------------------------------------------------------------
    // perform Cholesky factorization for A
    //--------------------------------------------------------------------------

    clock_t start = clock();

    DEMO_OK(SPEX_cholesky_analyze(&S, A, option));
    DEMO_OK(SPEX_cholesky_factorize(&F, A, S, option));

    clock_t end = clock();

    double t= (double) (end - start) / CLOCKS_PER_SEC;

    printf("\nSPEX Cholesky Factor time: %lf\n", t);

    //--------------------------------------------------------------------------
    // convert factorization to be updatable
    //--------------------------------------------------------------------------

    start = clock();

    DEMO_OK(SPEX_factorization_convert(F, true, option));

    end = clock();

    t= (double) (end - start) / CLOCKS_PER_SEC;

    printf("\ntime to make factorization updatable: %lf\n", t);

    //--------------------------------------------------------------------------
    // create a n-by-1 SPEX_DYNAMIC_CSC SPEX_MPZ matrix with w(0,0) = 1,
    // which will be used to update A as A = A + sigma*w*w^T
    //--------------------------------------------------------------------------
    // allocate an empty n-by-1 SPEX_DYNAMIC_CSC MPZ matrix
    DEMO_OK(SPEX_matrix_allocate(&w, SPEX_DYNAMIC_CSC, SPEX_MPZ, A->m, 1, 0, false,
        true, option));
    // reallocate w->v[0] with 1 entry
    DEMO_OK(SPEX_vector_realloc(w->v[0], 1, option));
    // set w(0,0) = 1
    DEMO_OK(SPEX_mpz_set_ui(w->v[0]->x[0], 1));
    w->v[0]->i[0] = 0;
    w->v[0]->nz = 1;

    //--------------------------------------------------------------------------
    // perform rank-1 Cholesky update/downdate
    //--------------------------------------------------------------------------

    start = clock();
    // compute the factorization of A = A + sigma*w*w^T by updating F
    int64_t sigma = 1; // sigma > 0 for update and sigma < 0 for downdate
    DEMO_OK(SPEX_update_cholesky_rank1(F, w, sigma, option));

    end = clock();

    t= (double) (end - start) / CLOCKS_PER_SEC;

    printf("\ntime to update factorization: %lf\n", t);

    //--------------------------------------------------------------------------
    // free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;
    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}
