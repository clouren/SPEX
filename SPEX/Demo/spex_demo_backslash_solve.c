//------------------------------------------------------------------------------
// Demo/spex_demo_backslash: example of SPEX_Solve for different factorizations
//------------------------------------------------------------------------------

// SPEX: (c) 2021-2025, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// A demo of SPEX_solve in C
 
# include "spex_demos.h"

#define FREE_WORKSPACE                           \
{                                                \
    if (mat_file != NULL)                        \
    {                                            \
        fclose(mat_file);                        \
    }                                            \
    mat_file = NULL ;                            \
    if (rhs_file != NULL)                        \
    {                                            \
        fclose(rhs_file);                        \
    }                                            \
    rhs_file = NULL ;                            \
    SPEX_matrix_free(&A,NULL);                   \
    SPEX_matrix_free(&AT,NULL);                  \
    SPEX_matrix_free(&b,NULL);                   \
    SPEX_matrix_free(&x_LU,NULL);                \
    SPEX_matrix_free(&x_CHOL,NULL);              \
    SPEX_matrix_free(&x_LDL,NULL);               \
    SPEX_matrix_free(&x_LU_T,NULL);              \
    SPEX_matrix_free(&x_CHOL_T,NULL);            \
    SPEX_matrix_free(&x_LDL_T,NULL);             \
    SPEX_factorization_free(&F_LU,NULL);         \
    SPEX_factorization_free(&F_CHOL,NULL);       \
    SPEX_factorization_free(&F_LDL,NULL);        \
    SPEX_symbolic_analysis_free(&S_LU, option);  \
    SPEX_symbolic_analysis_free(&S_LDL, option); \
    SPEX_symbolic_analysis_free(&S_CHOL, option);\
    SPEX_FREE(option);                           \
    SPEX_finalize();                             \
}                                                \

int main( int argc, char *argv[] )
{
    int64_t n = 0 ;
    SPEX_matrix A = NULL;
    SPEX_matrix AT = NULL;
    SPEX_matrix b = NULL;
    FILE *mat_file = NULL ;
    FILE *rhs_file = NULL ;
    SPEX_symbolic_analysis S_CHOL = NULL;
    SPEX_symbolic_analysis S_LDL = NULL;
    SPEX_symbolic_analysis S_LU = NULL;
    SPEX_factorization F_LU = NULL;
    SPEX_factorization F_LDL = NULL;
    SPEX_factorization F_CHOL = NULL;
    SPEX_matrix x_LU = NULL;
    SPEX_matrix x_LDL = NULL;
    SPEX_matrix x_CHOL = NULL;
    SPEX_matrix x_LU_T = NULL;
    SPEX_matrix x_LDL_T = NULL;
    SPEX_matrix x_CHOL_T = NULL;
    SPEX_options option = NULL;
    char *mat_name = NULL, *rhs_name = NULL;
    int64_t rat = 1;

    //--------------------------------------------------------------------------
    // Prior to using SPEX, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_TRY (SPEX_initialize ( )) ;

    // Set default options
    SPEX_TRY (SPEX_create_default_options(&option));

    // Process the command line
    SPEX_TRY (spex_demo_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    // Read in A
    mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return (1) ;
    }

    SPEX_TRY (spex_demo_tripread(&A, mat_file, SPEX_FP64, option));
    fclose(mat_file);
    mat_file = NULL ;
    n = A->n;
    // For this code, we utilize a vector of all ones as the RHS vector
    SPEX_TRY (SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option));

    // Read in b. The output of this demo function is b in dense format with
    // mpz_t entries
    rhs_file = fopen(rhs_name,"r");
    if( rhs_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return (1) ;
    }
    SPEX_TRY (spex_demo_read_dense(&b, rhs_file, option));
    fclose(rhs_file);
    rhs_file = NULL ;

    // Create A transpose in order to check the transpose solves
    SPEX_TRY( SPEX_transpose(&AT, A, option));

    //--------------------------------------------------------------------------
    // The demo reads in a SPD matrix. We will attempt an LU, LDL, and Cholesky
    // factorization and the general purpose solve
    //--------------------------------------------------------------------------

    option->order = SPEX_COLAMD;
    SPEX_TRY(SPEX_lu_analyze(&S_LU, A, option));

    option->order = SPEX_AMD;
    SPEX_TRY(SPEX_cholesky_analyze(&S_CHOL, A, option));
    SPEX_TRY(SPEX_ldl_analyze(&S_LDL, A, option));

    printf("checking SPEX_solve with LU ...\n");
    SPEX_TRY (SPEX_lu_factorize(&F_LU, A, S_LU, option));

    SPEX_TRY (SPEX_solve(&x_LU, F_LU, b, option));

    option->print_level=1;
    SPEX_TRY ( spex_demo_check_solution(A,x_LU,b,option));

    printf("checking SPEX_solve with Chol ...\n");

    SPEX_TRY (SPEX_cholesky_factorize(&F_CHOL, A, S_CHOL, option));

    SPEX_TRY (SPEX_solve(&x_CHOL, F_CHOL, b, option));

    SPEX_TRY ( spex_demo_check_solution(A,x_CHOL,b,option));

    printf("checking SPEX_solve with LDL ...\n");

    SPEX_TRY (SPEX_ldl_factorize(&F_LDL, A, S_LDL, option));

    SPEX_TRY (SPEX_solve(&x_LDL, F_LDL, b, option));

    SPEX_TRY ( spex_demo_check_solution(A,x_LDL,b,option));

    printf("checking SPEX_transpose_solve with LU ...\n");
    SPEX_TRY (SPEX_tsolve(&x_LU_T, F_LU, b, option));

    SPEX_TRY (spex_demo_check_solution(AT, x_LU_T, b, option));

    printf("checking SPEX_transpose solve with Cholesky ...\n");
    SPEX_TRY (SPEX_tsolve(&x_CHOL_T, F_CHOL, b, option));

    SPEX_TRY (spex_demo_check_solution(AT, x_CHOL_T, b, option));

    printf("checking SPEX_transpose solve with LDL ...\n");
    SPEX_TRY (SPEX_tsolve(&x_LDL_T, F_LDL, b, option));

    SPEX_TRY (spex_demo_check_solution(AT, x_CHOL_T, b, option));

    printf("\nAll SPEX_Solve tests successful!\n");

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
    return (0) ;
}

