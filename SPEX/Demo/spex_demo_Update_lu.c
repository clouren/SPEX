//------------------------------------------------------------------------------
// Demo/spex_update_demo_lu.c: demo for SPEX_Update library
//------------------------------------------------------------------------------

// SPEX: (c) 2020-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* A simple example to show how to perform LU factorization update for column
 * replacement.
 */

#include "spex_demos.h"

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

    SPEX_initialize ();

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    SPEX_options option = NULL;
    SPEX_factorization F = NULL;
    SPEX_symbolic_analysis S = NULL;
    SPEX_matrix A = NULL, A_DCSC = NULL;
    SPEX_matrix vk = NULL;
    SPEX_TRY (SPEX_create_default_options(&option));

    //--------------------------------------------------------------------------
    // read matrix and store as a SPEX_CSC SPEX_MPZ matrix A
    //--------------------------------------------------------------------------

    char *mat_name, *rhs_name;
    int64_t rat = 1;
    SPEX_TRY (spex_demo_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));
    
    FILE *mat_file = fopen(mat_name, "r");
    if (mat_file == NULL)
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    SPEX_TRY (spex_demo_tripread(&A, mat_file, SPEX_MPZ ,option));
    fclose(mat_file);

    //--------------------------------------------------------------------------
    // perform LU factorization for A
    //--------------------------------------------------------------------------

    double start = SUITESPARSE_TIME ;

    SPEX_TRY (SPEX_lu_analyze(&S, A, option));
    SPEX_TRY (SPEX_lu_factorize(&F, A, S, option));

    double end = SUITESPARSE_TIME ;

    double t= (end - start) ;

    printf("\nSPEX Left LU Factor time: %lf\n", t);

    //--------------------------------------------------------------------------
    // convert factorization to be updatable
    //--------------------------------------------------------------------------

    start = SUITESPARSE_TIME ;

    SPEX_TRY (SPEX_factorization_convert(F, true, option));

    end = SUITESPARSE_TIME ;

    t= (end - start) ;

    printf("\ntime to make factorization updatable: %lf\n", t);

    //--------------------------------------------------------------------------
    // create a n-by-1 SPEX_DYNAMIC_CSC SPEX_MPZ matrix with vk(0,0) = 1,
    // which will be used to replace one column from A
    //--------------------------------------------------------------------------
    // allocate an empty n-by-1 SPEX_DYNAMIC_CSC MPZ matrix
    SPEX_TRY (SPEX_matrix_allocate(&vk, SPEX_DYNAMIC_CSC, SPEX_MPZ, A->m, 1, 0, false,
        true, option));
    // reallocate vk->v[0] with 1 entry
    SPEX_TRY (SPEX_vector_realloc(vk->v[0], 1, option));
    // set vk(0,0) = 1
    SPEX_TRY (SPEX_mpz_set_ui(vk->v[0]->x[0], 1));
    vk->v[0]->i[0] = 0;
    vk->v[0]->nz = 1;

    //--------------------------------------------------------------------------
    // perform LU update for column replacement
    //--------------------------------------------------------------------------

    start = SUITESPARSE_TIME ;
    // compute the factorization of A after replacing the k-th column of A
    // with vk->v[0]
    int64_t k = 0;
    SPEX_TRY (SPEX_update_lu_colrep(F, vk, k, option));

    end = SUITESPARSE_TIME ;

    t= (end - start) ;

    printf("\ntime to update factorization: %lf\n", t);

    //--------------------------------------------------------------------------
    // optional: obtain updated matrix A
    //--------------------------------------------------------------------------

    // obtain A in SPEX_DYNAMIC_CSC MPZ
    SPEX_TRY (SPEX_matrix_copy(&A_DCSC, SPEX_DYNAMIC_CSC, SPEX_MPZ, A, option));
    // get the updated matrix A
    SPEX_TRY (SPEX_update_matrix_colrep(A_DCSC, vk, 0, option));

    //--------------------------------------------------------------------------
    // free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;
    printf ("\n%s: all tests passed\n\n", __FILE__);
    return 0;
}

