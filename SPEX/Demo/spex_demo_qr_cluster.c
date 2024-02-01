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
    int64_t n = 0, ok,m;

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
    //option->print_level = 4;
    DEMO_OK(spex_demo_tripread(&A, mat_file, SPEX_FP64, option));
    fclose(mat_file);

    n = A->n;
    m= A->m;
    // For this code, we utilize a vector of all ones as the RHS vector
    SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, m, 1, 0, false, true, option);
    if(!b){
        printf("scream \n");
    }

    // Create RHS
    for (int64_t k = 0; k < m; k++)
        DEMO_OK(SPEX_mpz_set_ui(b->x.mpz[k],1));


    //--------------------------------------------------------------------------
    // Dense test
    //--------------------------------------------------------------------------
  /*  SPEX_matrix A2 = NULL;     // Matrix to be randomly generated
    SPEX_matrix Q2 = NULL;
    SPEX_matrix R2 = NULL;
    SPEX_matrix b_new = NULL;
    SPEX_matrix x2 = NULL;

    DEMO_OK (SPEX_qr_analyze(&S, A, option));

    SPEX_matrix_copy(&A2, SPEX_DENSE, SPEX_MPZ, A, option);
    clock_t start_f = clock();
    SPEX_QR_IPGE( A2, &R2, &Q2);
    clock_t end_f= clock();

    option->print_level = 3;
    SPEX_matrix_check(Q2, option);

    double t_s = (double) (start_f)/CLOCKS_PER_SEC;
    double t_e= (double) (end_f)/CLOCKS_PER_SEC;
    printf("%s  %f %f\n", mat_name,  t_s, t_e);
    /*clock_t start_s = clock();
    SPEX_Qtb(Q2, b, &b_new);

    SPEX_QR_backsolve(R2, b_new, &x2);
    clock_t end_s = clock();

    double t_f = (double) (end_f - start_f) / CLOCKS_PER_SEC;
    double t_s =  (double) (end_s - start_s) / CLOCKS_PER_SEC;
    printf("%s %ld %ld  %f %f\n", mat_name, A->n, A->m, t_f, t_s);*/

    //--------------------------------------------------------------------------
    // Perform Analysis of A
    //--------------------------------------------------------------------------
    clock_t start_col = clock();
    //option->order =  SPEX_NO_ORDERING;
    DEMO_OK (SPEX_qr_analyze(&S, A, option));
    clock_t end_col = clock();

    //--------------------------------------------------------------------------
    // Factorize AQ
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();
    option->print_level = 3;
    DEMO_OK (SPEX_qr_factorize(&F, A, S, option));
    clock_t end_factor = clock();
    //SPEX_matrix_check(F->Q, option);
    //SPEX_matrix_check(F->R, option);


    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();
    DEMO_OK (SPEX_qr_solve(&x, F, b, option));
    clock_t end_solve = clock();
    //SPEX_matrix_check(x, option);
     
    double t_sym = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;
    printf("%s %ld %ld %ld %ld %ld %ld %f %f %f\n", mat_name, A->n, A->m, F->rank, A->p[A->n], (F->Q->p[F->Q->n]), (F->R->p[F->R->n]), t_sym, t_factor, t_solve);

    // Check solution
    //option->print_level=1;
    //DEMO_OK(spex_demo_check_solution(A,x,b,option)); //works is x is mpq
   
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;                        
}

