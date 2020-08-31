//------------------------------------------------------------------------------
// SPEX_Chol/Demo/example.c: example main program for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M University.  All Rights Reserved.  
// See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------


/* This example shows how to use SPEX Chol with a given input matrix and a double
   output. The input is read from a file */

// usage:
// example > out
// out is file for output calculated result

#define FREE_WORKSPACE              \
    SPEX_LU_analysis_free(&S, option);\
    SPEX_matrix_free(&A, option);   \
    SPEX_FREE(option);              \
    SPEX_matrix_free(&b, option);   \
    SPEX_matrix_free(&x, option);   \
    SPEX_finalize();

#include "SPEX_Chol.h"   
    
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
    mat_name = "../ExampleMats/2.mat";
    rhs_name = "../ExampleMats/2.mat.soln";
    if (argc > 2)
    {
        mat_name = argv[1];
        rhs_name = argv[2];
    }

    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_info ok;
    SPEX_matrix *A = NULL ;                     // input matrix
    SPEX_matrix *b = NULL ;                     // Right hand side vector
    SPEX_matrix *x = NULL ;                     // Solution vectors
    SPEX_LU_analysis *S = NULL ;                // Column permutation
    SPEX_options *option = SPEX_create_default_options();
    if (option == NULL)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Allocate memory, read in A and b
    //--------------------------------------------------------------------------

    // Read in A. The output of this demo function is A in CSC format with
    // mpz_t entries.
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    OK(SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);

    // Read in b. The output of this demo function is b in dense format with
    // mpz_t entries
    //FILE* rhs_file = fopen(rhs_name,"r");
    //if( rhs_file == NULL )
    //{
    //    perror("Error while opening the file");
    //    FREE_WORKSPACE;
    //    return 0;
    //}
//    OK(SPEX_read_dense(&b, rhs_file, option));
    int64_t n = A->n;
    SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true, option);
    // Create RHS
    for (int64_t k = 0; k < n; k++)
        OK(SPEX_mpz_set_ui(b->x.mpz[k],1));

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
    
    // SPEX LU has an optional check, to enable it, one can set the following
    // parameter to be true.
    option->check = true;
   
    // Solve the system and give MPQ solution
    OK(SPEX_Chol_backslash( &x, SPEX_MPQ, A, b, option));
    
    clock_t end_s = clock();

    double t_s = (double) (end_s - start_s) / CLOCKS_PER_SEC;

    printf("\nSPEX Chol Factor & Solve time: %lf\n", t_s);

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;

    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}

