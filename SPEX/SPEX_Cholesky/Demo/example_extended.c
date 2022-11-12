//------------------------------------------------------------------------------
// SPEX_Cholesky/Demo/example_extended: Demo main program for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy,
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Include the Integer-preserving Cholesky routines */

#define FREE_WORKSPACE                          \
{                                               \
    SPEX_matrix_free(&A,NULL);                  \
    SPEX_matrix_free(&b,NULL);                  \
    SPEX_matrix_free(&x,NULL);                  \
    SPEX_symbolic_analysis_free (&S, option);   \
    SPEX_matrix_free(&((F)->L), option);        \
    SPEX_factorization_free(&F, option);        \
    SPEX_FREE(option);                          \
    SPEX_finalize();                            \
}

#define DEMO_OK(method)                         \
{                                               \
    ok = method ;                               \
    if (ok != SPEX_OK)                          \
    {                                           \
        SPEX_cholesky_determine_error(ok);      \
        FREE_WORKSPACE ;                        \
        return 0 ;                              \
    }                                           \
}

#include "chol_demos.h"

int main( int argc, char *argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX-Chol, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int64_t n = 0, ok;

    SPEX_symbolic_analysis S = NULL;
    SPEX_factorization F = NULL ;
    SPEX_matrix A = NULL;
    SPEX_matrix b = NULL;
    SPEX_matrix x = NULL;

    // Default options.
    SPEX_options option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));

    // FIXME: Set demo matrix and RHS name
    char* mat_name = "../../ExampleMats/2.mat.txt";
    char* rhs_name = "../../ExampleMats/2.mat.soln.txt";
    int64_t rat = 1;

    // Process the command line
    DEMO_OK(SPEX_cholesky_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    // Read in A
    FILE *mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }

    DEMO_OK(SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);
    n = A->n;
    // For this code, we utilize a vector of all ones as the RHS vector
    SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option);
    // Create RHS
    for (int64_t k = 0; k < n; k++)
        DEMO_OK(SPEX_mpz_set_ui(b->x.mpz[k],1));

    //--------------------------------------------------------------------------
    // Perform Analysis of A
    //--------------------------------------------------------------------------

    clock_t start_col = clock();

    // Symmetric ordering of A. Uncomment the desired one, AMD is recommended
    //option->order = SPEX_NO_ORDERING;  // No ordering
    option->order = SPEX_AMD;  // AMD
    //option->order = SPEX_COLAMD; // COLAMD
    DEMO_OK(SPEX_cholesky_analyze(&S, A, option));
    //DEMO_OK(SPEX_cholesky_preorder(&S, A, option));
    clock_t end_col = clock();

    //--------------------------------------------------------------------------
    // Factorize PAP
    //--------------------------------------------------------------------------

    //option->algo=SPEX_CHOL_LEFT;
    clock_t start_factor = clock();

    DEMO_OK( SPEX_cholesky_factorize(&F, A, S, option));

    clock_t end_factor = clock();

    option->print_level=3;
    //DEMO_OK(SPEX_matrix_check(F->L,option));

    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------

    clock_t start_solve = clock();

    DEMO_OK( SPEX_cholesky_solve(&x, F, b, option));

    clock_t end_solve = clock();

    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------

    double t_col = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;

    printf("\nNumber of L nonzeros: \t\t\t%ld",
        (F->L->p[F->L->n]) );
    printf("\nSymbolic Analysis Check time: \t\t%lf", t_col);
    printf("\nIP Chol Factorization time: \t\t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);

    // Check solution
    option->print_level=1;
    DEMO_OK ( SPEX_check_solution(A,x,b,option));

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;
}

