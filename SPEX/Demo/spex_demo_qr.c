//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR_dense.c: Dense REF QR factorization
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

//(future) TODO:
/* 
* python + matlab interface
* user guide
* tcov
* tests with rank deficient (ideally also ill conditioned so we can compare to matlab?)
* tests with ls probs
* matrix multiply (to compare with chol for ls probs)
* actual demo like chol's
*/

# include "SPEX.h"
# include "spex_util_internal.h"
# include "spex_qr_internal.h"
# include "spex_cholesky_internal.h"
#include "spex_demos.h"

#define FREE_WORKSPACE                          \
    SPEX_matrix_free(&A,  NULL);                \
    SPEX_FREE(option);                          \
    SPEX_finalize(); 
    

    /*    SPEX_matrix_free(&A2, NULL);                \
    SPEX_matrix_free(&R2, NULL);                \
    SPEX_matrix_free(&Q2, NULL);                \*/
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

    // Input arguments.
    unsigned int seed;
    int64_t m;
    int64_t n;
    int64_t lower;
    int64_t upper;

    // Process the command line. If there are not exactly 5 input arguments,
    // then default values are used.
    if (argc != 6)
    {
        printf("\nExpected usage: ./SPEX_QR_dense SEED M N LOWER UPPER\n");
        printf("\nUsing default settings\n");
        seed = 10;
        //m = 100;
        m= 50;
        n = 50;
        lower = 1;
        upper = 10;
    }
    // Acquire input arguments
    else
    {
        seed = (unsigned int) atoi(argv[1]);
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        lower = atoi(argv[4]);
        upper = atoi(argv[5]);
    }

m=5;n=4;seed=14;
//colamd m<n
    // Input checks
    ASSERT(m >= 0);
    ASSERT(n >= 0);
    ASSERT(lower < upper);
    ASSERT(seed >= 1);

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    SPEX_info ok;
    SPEX_matrix A = NULL ;     // Integer matrix to be factorized
    SPEX_matrix A2 = NULL;     // Matrix to be randomly generated

    SPEX_matrix Ainit = NULL;     // Matrix to be randomly generated

    // Next we define 3 Q R pairs. Each pair is generated via a different dense
    // algorithm
    SPEX_matrix Q2 = NULL;
    SPEX_matrix R2 = NULL;

    SPEX_matrix b = NULL;
    SPEX_matrix b2 = NULL;
    SPEX_matrix b_new = NULL;
    SPEX_matrix x = NULL;
    SPEX_matrix x2 = NULL;
    
    // SPEX Options
    SPEX_options option = NULL;
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
/*
    SPEX_generate_random_matrix ( &Ainit, m, n, seed, lower, upper);
    Ainit->nz = m*n;
    option->print_level = 3;
    //SPEX_matrix_check(Ainit, option);

    // Create A as a copy of Ainit
    // A is a copy of the Ainit matrix. A is a sparse matrix with mpz_t entries
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, Ainit, option);

    // Create A2 as a copy of Ainit
    // A2 is a copy of the Ainit matrix. A is a dense matrix with mpz_t entries
    SPEX_matrix_copy(&A2, SPEX_DENSE, SPEX_MPZ, Ainit, option);
    
     option->print_level = 3;
     //SPEX_matrix_check(A, option);
     
     SPEX_generate_random_matrix ( &b2, m, 1, seed, lower, upper);
    b2->nz = m;
    // Make a copy of b
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b2, option);

    //option->print_level = 3;
    //SPEX_matrix_check(A, option);
/**/
    //char *mat_name = "ExampleMats/smallRankDeficient.mat.txt";
    //char *rhs_name = "ExampleMats/smallRankDeficient.rhs.txt";
    //char *mat_name = "ExampleMats/LF10.mat.txt";
    //char *rhs_name = "ExampleMats/LF10.rhs.txt";
    char *mat_name = "ExampleMats/smallZeros.mat.txt";
    char *rhs_name = "ExampleMats/smallZeros.rhs.txt";
    // Read in A
    FILE *mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the mat file");
        FREE_WORKSPACE;
        return 0;
    }

    DEMO_OK(spex_demo_tripread(&A, mat_file, SPEX_FP64, option));
    fclose(mat_file);
    n = A->n;
    m = A->m;

    // Read in b. The output of this demo function is b in dense format with
    // mpz_t entries
    FILE *rhs_file = fopen(rhs_name,"r");
    if( rhs_file == NULL )
    {
        perror("Error while opening the rhs file");
        FREE_WORKSPACE;
        return 0;
    }
    DEMO_OK(spex_demo_read_dense(&b, rhs_file, option));
    fclose(rhs_file);
    
    option->print_level = 3;
    //SPEX_matrix_check(A, option);*/

    //--------------------------------------------------------------------------
    // Dense
    //--------------------------------------------------------------------------
    /*
    SPEX_matrix_copy(&A2, SPEX_DENSE, SPEX_MPZ, A, option);
    option->print_level = 3;
    SPEX_QR_IPGE( A2, &R2, &Q2);
    //printf("Dense Q, R\n");
    //SPEX_matrix_check(Q2, option);
    //SPEX_matrix_check(R2, option);

    SPEX_Qtb(Q2, b, &b_new);
    //spex_matrix_mul(b_new,R->x.mpz[R->nz]);

    SPEX_QR_backsolve(R2, b_new, &x2);
    printf("orint x dense:\n");
    SPEX_matrix_check(x2, option);
    /*SPEX_mpq_set_num(x->scale, SPEX_2D(R, n-1, n-1, mpz));


    // Create a double version of x
    SPEX_matrix_copy(&x_doub, SPEX_DENSE, SPEX_FP64, x, option);
    SPEX_matrix_check(b_new, option);*/

    //--------------------------------------------------------------------------
    // Sparse
    //--------------------------------------------------------------------------

/**/
    SPEX_info info;
    SPEX_matrix rhos = NULL,R3=NULL, rhos2 = NULL;
    int64_t *h;
    int64_t j=0, nnz=n*m;
    int64_t i,pQ, pR;
    int sgn;
    SPEX_symbolic_analysis S = NULL;
    SPEX_factorization F = NULL ;

    printf("analysis:\n");
    //option->print_level = 3;
    option->order =  SPEX_NO_ORDERING;
    DEMO_OK (SPEX_qr_analyze(&S, A, option));
    SPEX_matrix_check(A, option); 
    printf("facts:\n");
    option->print_level = 3;
    DEMO_OK (SPEX_qr_factorize(&F, A, S, option));
    SPEX_matrix_check(F->Q, option);
    SPEX_matrix_check(F->R, option);

    printf("solve:\n");
    DEMO_OK (SPEX_qr_solve(&x, F, b, option));
    printf("Success!!\n");

    /*
    //option->order =  SPEX_NO_ORDERING;
    //SPEX_qr_backslash(&x,SPEX_FP64,A,b, option);
    */

    //printf("orint x sparse:\n");
     SPEX_matrix_check(x, option);
     
     printf("Rank of matrix: %ld, is deficient? %ld\n",F->rank,(F->R->n)-(F->rank)); //if rank defficient then sol will be wrong
     DEMO_OK(spex_demo_check_solution(A,x,b,option)); //works is x is mpq
   /* */
    ////
    // Tests
    ////
    /*
    SPEX_info info;
    SPEX_info ok;
    
    char *mat_name, *rhs_name;
    int64_t rat = 1;
    
    SPEX_matrix A = NULL ; 
    SPEX_symbolic_analysis S = NULL;
    SPEX_factorization F = NULL ;
    SPEX_options option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));
    if (option == NULL)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    
    //DEMO_OK(spex_demo_process_command_line(argc, argv, option,
       // &mat_name, &rhs_name, &rat));
    mat_name = argv[2];
    //printf("%s\n",mat_name);

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
    
    DEMO_OK(spex_demo_tripread(&A, mat_file, SPEX_FP64, option));
    fclose(mat_file);
    
    int64_t n = A->n, col_sum;
    int sgn;
    option->order =  SPEX_NO_ORDERING;
    DEMO_OK (SPEX_qr_analyze(&S, A, option));

    DEMO_OK (SPEX_qr_factorize(&F, A, S, option));
    
    
    printf("%s, ",mat_name);
    printf("%ld,  ",n);
    printf("%ld, ",F->R->nz);
    
   /* for(int64_t i=0;i<n;i++)
    {
        col_sum=0;
        for(int64_t p=F->Q->p[i]; p<F->Q->p[i+1]; p++)
        {
            SPEX_MPZ_SGN(&sgn, F->Q->x.mpz[p]); //Q(i,k)
            if(sgn!=0)
            {
                col_sum++;
            }
            //printf("%ld %ld\n",i,p);
        }
        printf(" %ld, ", col_sum);
    }
    printf("\n");*/
/*
    printf("%ld,\n ",F->Q->nz);
    */
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    //FREE_WORKSPACE;
    SPEX_factorization_free(&F, NULL);      
    SPEX_symbolic_analysis_free(&S, NULL);    
    SPEX_matrix_free(&A,  NULL);              
    //SPEX_matrix_free(&x,  NULL);   
    SPEX_FREE(option);                         
    

}

