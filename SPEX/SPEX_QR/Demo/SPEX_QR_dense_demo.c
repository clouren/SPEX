//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR_dense.c: Dense REF QR factorization
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code contains a dense REF QR factorization. It is meant to be proof of 
 * concept for the REF QR factorization algorithms
 */


# include "SPEX_QR.h"

#define FREE_WORKSPACE                          \
    SPEX_matrix_free(&A,  NULL);                \
    SPEX_matrix_free(&A2, NULL);                \
    SPEX_matrix_free(&R,  NULL);                \
    SPEX_matrix_free(&Q,  NULL);                \
    SPEX_matrix_free(&R2, NULL);                \
    SPEX_matrix_free(&Q2, NULL);                \
    SPEX_matrix_free(&R3, NULL);                \
    SPEX_matrix_free(&Q3, NULL);                \
    SPEX_matrix_free(&b, NULL);                 \
    SPEX_matrix_free(&b2, NULL);                \
    SPEX_matrix_free(&b_new, NULL);             \
    SPEX_FREE(option);                          \
    SPEX_finalize();                            \


int main( int argc, char* argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX QR, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_initialize();
    
    unsigned int seed;
    int64_t m;
    int64_t n;
    int64_t lower;
    int64_t upper;
    
    if (argc != 6)
    {
        printf("\nExpected usage: ./SPEX_QR_dense SEED M N LOWER UPPER\n");
        printf("\nUsing default settings\n");
        seed = 10;
        m = 100;
        n = 100;
        lower = 1;
        upper = 10;
    }
    else
    {
        seed = (unsigned int) atoi(argv[1]);
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        lower = atoi(argv[4]);
        upper = atoi(argv[5]);
    }
    

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    SPEX_info ok;
    SPEX_matrix *A = NULL ;                     // input matrix
    SPEX_matrix *A2 = NULL;                    // Matrix to be generated
    SPEX_matrix *R = NULL;                      // Upper triangular matrix
    SPEX_matrix *Q = NULL;                      // Orthogonal Matrix
    SPEX_matrix *Q2 = NULL;
    SPEX_matrix *R2 = NULL;
    SPEX_matrix *Q3 = NULL;
    SPEX_matrix *R3 = NULL;
    SPEX_matrix *b = NULL;
    SPEX_matrix *b2 = NULL;
    SPEX_matrix *b_new = NULL;
    SPEX_options *option = NULL;
    SPEX_create_default_options(&option);
    if (!option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Generate a random dense matrix
    //--------------------------------------------------------------------------
    
    SPEX_generate_random_matrix ( &A2, m, n, seed, lower, upper);
    A2->nz = m*n;
    
    // Create A as a copy of A2
    // A is a copy of the A2 matrix. A is a dense matrix with mpz_t entries
    SPEX_matrix_copy(&A, SPEX_DENSE, SPEX_MPZ, A2, option);
    
    //--------------------------------------------------------------------------
    // Factorize
    //--------------------------------------------------------------------------
    
    // Perform the mulitplication heavy QR_IPGE. This method is focused on doing
    // dot products and tries to limit the number of divisions
    // Better for parallelization and memory
    // I.e., algorithm 1 from paper
    
    clock_t start_solve1 = clock();
    
    SPEX_QR_IPGE( A, &R, &Q);
    
    clock_t end_solve1 = clock();
    
    // Use IPGE Pursell method. This approach performs IPGE on [A' * A | A']
    // This specific version operates on A'*A and A' seperately
    
    
    clock_t start_solve2 = clock();
    
    SPEX_QR_PURSELL( A, &R2, &Q2);
    
    clock_t end_solve2 = clock();
    
    // Use IPGE Pursell method. This approach explicitly constructs [A'*A | A']
    // The pursell method is generally quite division heavy.
    
       
    clock_t start_solve3 = clock();
    
    SPEX_QR_PURSELL2( A, &R3, &Q3);
    
    clock_t end_solve3 = clock();

    //option->print_level = 3;
    //SPEX_matrix_check(Q, option);
    
    //SPEX_matrix_check(Q2, option);
    
    //SPEX_matrix_check(Q3, option);
    
    // Now check to make sure they are the same
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->m; j++)
        {
            int r ;
            // Have to transpose Q here because theirs is backwards
            SPEX_mpz_cmp(&r, SPEX_2D(Q,j,i,mpz), SPEX_2D(Q2, i, j, mpz));
            if ( r != 0)
                printf("\n Q2 Incorrect at %ld %ld", i, j);
        }
    }
    
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            SPEX_mpz_cmp(&r, SPEX_2D(R,i,j,mpz), SPEX_2D(R2, i, j, mpz));
            if ( r != 0)
                printf("\n R2 Incorrect at %ld %ld", i, j);
        }
    }
    
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            SPEX_mpz_cmp(&r, SPEX_2D(Q,i,j,mpz), SPEX_2D(Q3, i, j, mpz));
            if ( r != 0)
                printf("\n Q3 Incorrect at %ld %ld", i, j);
        }
    }
    
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            SPEX_mpz_cmp(&r, SPEX_2D(R,i,j,mpz), SPEX_2D(R3, i, j, mpz));
            if ( r != 0)
                printf("\n R3 Incorrect at %ld %ld", i, j);
        }
    }
    
    
    //--------------------------------------------------------------------------
    // Randomly generate a RHS
    //--------------------------------------------------------------------------
    
    // Generate floating point matrix
    SPEX_generate_random_matrix ( &b2, m, 1, seed, lower, upper);
    b2->nz = m;
    
    clock_t start_b = clock();
    
    // Make a copy of b
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b2, option);
    
    // Compute Q'*b
    SPEX_Qtb(Q, b, &b_new);
    
    clock_t end_b = clock();
    
    // TODO Solve x = R \ Q^T b
    
    
    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------
    
    double t_qr1 =  (double) (end_solve1 - start_solve1) / CLOCKS_PER_SEC;
    double t_qr2 =  (double) (end_solve2 - start_solve2) / CLOCKS_PER_SEC;
    double t_qr3 =  (double) (end_solve3 - start_solve3) / CLOCKS_PER_SEC;
    double t_b =  (double) (end_b - start_b) / CLOCKS_PER_SEC;


    double rat  = t_qr2 / t_qr1;
    double rat2 = t_qr3 / t_qr1;
    
    printf("\nIPGE QR time: \t\t\t%lf\n\n", t_qr1);
    printf("\nIPGE QR Pursell time: \t\t%lf\n\n", t_qr2);
    printf("\nIPGE QR Pursell2 time: \t\t%lf\n\n", t_qr3);
    printf("\nRatio IPGE QR/ Pursell1: \t%lf\n\n", rat);
    printf("\nRatio IPGE QR/ Pursell2: \t%lf\n\n", rat2);
    printf("\nTime to do Q^T*b: \t\t%lf\n", t_b);
        
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
                 
}
    
