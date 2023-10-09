//------------------------------------------------------------------------------
// SPEX_Demo/spex_demo_qr_extended.c: 
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_demos.h"

#define FREE_WORKSPACE                          \
{                                               \
    SPEX_matrix_free(&A,NULL);                  \
    SPEX_matrix_free(&b,NULL);                  \
    SPEX_matrix_free(&x,NULL);                  \
    SPEX_symbolic_analysis_free (&S, option);   \
    SPEX_factorization_free(&F, option);        \
    SPEX_FREE(option);                          \
    SPEX_finalize();                            \
}
    
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

    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int64_t n = 0, ok;

    SPEX_symbolic_analysis S = NULL;
    SPEX_factorization F = NULL ;
    SPEX_matrix A = NULL;
    SPEX_matrix b = NULL;
    SPEX_matrix x = NULL;
   
    char *mat_name, *rhs_name;
    int64_t rat = 1;

    // Default options.
    SPEX_options option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));

    // Process the command line
    DEMO_OK(spex_demo_process_command_line(argc, argv, option,
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

    DEMO_OK(spex_demo_tripread(&A, mat_file, SPEX_FP64, option));
    fclose(mat_file);
    n = A->n;
    // For this code, we utilize a vector of all ones as the RHS vector
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
    DEMO_OK(spex_demo_read_dense(&b, rhs_file, option));
    fclose(rhs_file);

    //--------------------------------------------------------------------------
    // Perform Analysis of A
    //--------------------------------------------------------------------------
    printf("Analysis:\n");
    DEMO_OK (SPEX_qr_analyze(&S, A, option));

    //--------------------------------------------------------------------------
    // Factorize AQ
    //--------------------------------------------------------------------------
    printf("Factorization:\n");
    DEMO_OK (SPEX_qr_factorize(&F, A, S, option));
    //SPEX_matrix_check(F->Q, option);
    //SPEX_matrix_check(F->R, option);

    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    printf("Solve:\n");
    DEMO_OK (SPEX_qr_solve(&x, F, b, option));
    //SPEX_matrix_check(x, option);
     
    printf("Success!!\n");
    printf("Rank of matrix: %ld, is deficient? %ld\n",F->rank,(F->R->n)-(F->rank));
    // Check solution
    option->print_level=1;
    DEMO_OK(spex_demo_check_solution(A,x,b,option)); //works is x is mpq
   
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;                        
}

