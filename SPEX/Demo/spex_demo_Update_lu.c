//------------------------------------------------------------------------------
// Demos/SPEX_Update_demo_lu.c: demo for SPEX_Update library
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* A simple example to show how to perform LU factorization update for column
 * replacement.
 */

#include "demos.h"

#define FREE_WORKSPACE                           \
{                                                \
    SPEX_matrix_free(&A, option);                \
    SPEX_matrix_free(&A_DCSC, option);           \
    SPEX_matrix_free(&vk, option);               \
    SPEX_symbolic_analysis_free(&S, option);     \
    SPEX_factorization_free(&F, option);         \
    SPEX_FREE(option);                           \
    SPEX_finalize();                             \
}

int main(int argc, char *argv[] )
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
    SPEX_matrix A = NULL, A_DCSC = NULL;
    SPEX_matrix vk = NULL;
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
    DEMO_OK(SPEX_tripread(&A, mat_file, option));
    fclose(mat_file);

    //--------------------------------------------------------------------------
    // perform LU factorization for A
    //--------------------------------------------------------------------------

    clock_t start = clock();

    DEMO_OK(SPEX_lu_analyze(&S, A, option));
    DEMO_OK(SPEX_lu_factorize(&F, A, S, option));

    clock_t end = clock();

    double t= (double) (end - start) / CLOCKS_PER_SEC;

    printf("\nSPEX Left LU Factor time: %lf\n", t);

    //--------------------------------------------------------------------------
    // convert factorization to be updatable
    //--------------------------------------------------------------------------

    start = clock();

    DEMO_OK(SPEX_factorization_convert(F, true, option));

    end = clock();

    t= (double) (end - start) / CLOCKS_PER_SEC;

    printf("\ntime to make factorization updatable: %lf\n", t);

    //--------------------------------------------------------------------------
    // create a n-by-1 SPEX_DYNAMIC_CSC SPEX_MPZ matrix with vk(0,0) = 1,
    // which will be used to replace one column from A
    //--------------------------------------------------------------------------
    // allocate an empty n-by-1 SPEX_DYNAMIC_CSC MPZ matrix
    DEMO_OK(SPEX_matrix_allocate(&vk, SPEX_DYNAMIC_CSC, SPEX_MPZ, A->m, 1, 0, false,
        true, option));
    // reallocate vk->v[0] with 1 entry
    DEMO_OK(SPEX_vector_realloc(vk->v[0], 1, option));
    // set vk(0,0) = 1
    DEMO_OK(SPEX_mpz_set_ui(vk->v[0]->x[0], 1));
    vk->v[0]->i[0] = 0;
    vk->v[0]->nz = 1;

    //--------------------------------------------------------------------------
    // perform LU update for column replacement
    //--------------------------------------------------------------------------

    start = clock();
    // compute the factorization of A after replacing the k-th column of A
    // with vk->v[0]
    int64_t k = 0;
    DEMO_OK(SPEX_update_lu_colrep(F, vk, k, option));

    end = clock();

    t= (double) (end - start) / CLOCKS_PER_SEC;

    printf("\ntime to update factorization: %lf\n", t);

    //--------------------------------------------------------------------------
    // optional: obtain updated matrix A
    //--------------------------------------------------------------------------

    // obtain A in SPEX_DYNAMIC_CSC MPZ
    DEMO_OK(SPEX_matrix_copy(&A_DCSC, SPEX_DYNAMIC_CSC, SPEX_MPZ, A, option));
    // get the updated matrix A
    DEMO_OK(SPEX_update_matrix_colrep(A_DCSC, vk, 0, option));

    //--------------------------------------------------------------------------
    // free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;
    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}

