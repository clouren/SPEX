//------------------------------------------------------------------------------
// SPEX_Chol/Demo/example.c: example main program for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

//TODO (maybe) change name of folder from SPEX_Cholesky to Cholesky
//TOASK valgrind: Conditional jump or move depends on uninitialised value(s) ??

/* This example shows how to use SPEX Chol with a given input matrix and a double
   output. The input is read from a file */

// usage:
// example > out
// out is file for output calculated result

#define FREE_WORKSPACE              \
    SPEX_matrix_free(&A, option);              \
    SPEX_matrix_free(&b, option);   \
    SPEX_matrix_free(&x, option);   \
    SPEX_matrix_free(&x2, option);   \
    SPEX_FREE(option);   \
    SPEX_finalize();

#include "demos.h"   

#define DEMO_OK(method)                 \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        SPEX_Chol_determine_error(ok);  \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}

    
int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // Prior to using SPEX Chol, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Get matrix and right hand side file names
    //--------------------------------------------------------------------------
    char *mat_name, *rhs_name;
    mat_name = "../ExampleMats/872.mat.txt";
    rhs_name = "../ExampleMats/872.mat.soln.txt";
    if (argc > 2)
    {
        mat_name = argv[1];
        rhs_name = argv[2];
    }

    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_info ok;
    //SPEX_matrix *Ad = NULL ;                     // input matrix as double
    SPEX_matrix* A = NULL ;                     // input matrix with mpz values
    SPEX_matrix* b = NULL ;                     // Right hand side vector
    SPEX_matrix* x = NULL ;                     // Solution vectors
    SPEX_matrix* x2 = NULL ;                     // copy of solution vectors
    SPEX_options* option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));
    if (option == NULL)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    option->order = SPEX_AMD; //AMD is default for Cholesky //TODO maybe change Util?
    //--------------------------------------------------------------------------
    // Allocate memory, read in A and b
    //--------------------------------------------------------------------------

    // Read in A. The output of this demo function is A in CSC format with
    // double entries.
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    
    DEMO_OK(SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);

    int64_t n = A->n;
    SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true, option);

    // Create RHS
    for (int64_t k = 0; k < n; k++)
        DEMO_OK(SPEX_mpz_set_ui(b->x.mpz[k],1));

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
    
    // SPEX Cholesky has an optional check, to enable it, one can set the following
    // parameter to be true.
    option->check = true;
    //option->print_level = 1;
    // Solve the system and give MPQ solution
    DEMO_OK(SPEX_Chol_backslash( &x, SPEX_MPQ, A, b, option));
    
    clock_t end_s = clock();

    double t_s = (double) (end_s - start_s) / CLOCKS_PER_SEC;

    printf("\nSPEX Chol Factor & Solve time: %lf\n", t_s);

    // x2 is a copy of the solution. x2 is a dense matrix with mpfr entries
    DEMO_OK ( SPEX_matrix_copy(&x2, SPEX_DENSE, SPEX_FP64, x, option));
    
    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;

    printf ("\n%s: all tests passed\n\n", __FILE__) ;

    return 0;
}

