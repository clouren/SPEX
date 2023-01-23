//------------------------------------------------------------------------------
// Demos/SPEX_Chol_demo_simple.c: example of simple call of SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* This example shows how to use SPEX Chol with a given input matrix and a
 * double output. The input is read from a file */

// usage:
// example > out
// out is file for output calculated result

#include "demos.h"

#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free(&A, option);           \
    SPEX_matrix_free(&b, option);           \
    SPEX_matrix_free(&x, option);           \
    SPEX_matrix_free(&x2, option);          \
    SPEX_FREE(option);                      \
    SPEX_finalize();                        \
}


int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // Prior to using SPEX Chol, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Get matrix file name
    //--------------------------------------------------------------------------

    char *mat_name;
    char *rhs_name; //this is actually ignored and we're using a rhs of 1s
    int64_t rat = 1;
    if (argc > 2)
    {
        mat_name = argv[1];
    }

    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_info ok;
    SPEX_matrix A = NULL ;                     // input matrix with mpz values
    SPEX_matrix b = NULL ;                     // Right hand side vector
    SPEX_matrix x = NULL ;                     // Solution vectors
    SPEX_matrix x2 = NULL ;                     // copy of solution vectors
    SPEX_options option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));
    if (option == NULL)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    option->order = SPEX_AMD; //AMD is default for Cholesky
    
    // Process the command line
    DEMO_OK(SPEX_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));

    //--------------------------------------------------------------------------
    // Allocate memory, read in A and b
    //--------------------------------------------------------------------------

    // Read in A. The output of this demo function is A in CSC format with
    // double entries.
    FILE *mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }

    DEMO_OK(SPEX_tripread(&A, mat_file, SPEX_FP64, option));
    fclose(mat_file);

    int64_t n = A->n;
    SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option);

    // Read in b. The output of this demo function is b in dense format with
    // mpz_t entries
    FILE *rhs_file = fopen(rhs_name,"r");
    if( rhs_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    DEMO_OK(SPEX_read_dense(&b, rhs_file, option));
    fclose(rhs_file);

    // Check if the size of A matches b
    if (A->n != b->m)
    {
        printf("%"PRId64" %"PRId64" \n", A->m,b->m);
        fprintf (stderr, "Error! Size of A and b do not match!\n");
        FREE_WORKSPACE;
        return 0;
    }
    //--------------------------------------------------------------------------
    // solve
    //--------------------------------------------------------------------------
    clock_t start_s = clock();
    option->algo=SPEX_CHOL_LEFT;

    DEMO_OK(SPEX_cholesky_backslash( &x, SPEX_MPQ, A, b, option));

    clock_t end_s = clock();

    double t_s = (double) (end_s - start_s) / CLOCKS_PER_SEC;

    printf("\nSPEX Chol Factor & Solve time: %lf\n", t_s);

    // x2 is a copy of the solution. x2 is a dense matrix with mpfr entries
    DEMO_OK ( SPEX_matrix_copy(&x2, SPEX_DENSE, SPEX_FP64, x, option));

    // check solution
    option->print_level=1;
    DEMO_OK ( SPEX_check_solution(A,x,b,option));

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;

    printf ("\n%s: all tests passed\n\n", __FILE__) ;

    return 0;
}

