//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR_dense.c: Dense REF QR factorization
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* This code contains a dense REF QR factorization. It is meant to be proof of
 * concept for the REF QR factorization algorithms. It computes the
 * factorization A = Q D R where Q and R is integer. Note that all code in this
 * version of SPEX QR assumes that A is a fully dense matrix; thus these
 * routines are not yet appropriate for sparse matrices.
 */


# include "SPEX.h"
# include "spex_util_internal.h"
# include "spex_qr_internal.h"
# include "spex_cholesky_internal.h"
#include "spex_demos.h"

#define FREE_WORKSPACE                          \
    SPEX_matrix_free(&A,  NULL);                \
    SPEX_matrix_free(&A2, NULL);                \
    SPEX_matrix_free(&R2, NULL);                \
    SPEX_matrix_free(&Q2, NULL);                \
    SPEX_FREE(option);                          \
    SPEX_finalize();                            \

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

m=5;n=5;seed=14;
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

    /*SPEX_generate_random_matrix ( &Ainit, m, n, seed, lower, upper);
    Ainit->nz = m*n;

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
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b2, option);*/
    char *mat_name = "ExampleMats/smallZeros.mat.txt";
    char *rhs_name = "ExampleMats/smallZeros.rhs.txt";
    //char *mat_name = "ExampleMats/LF10.mat.txt";
    //char *rhs_name = "ExampleMats/LF10.rhs.txt";
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
    
    //option->print_level = 3;
    //SPEX_matrix_check(A, option);

    //--------------------------------------------------------------------------
    // Dense
    //--------------------------------------------------------------------------
    
   /*SPEX_matrix_copy(&A2, SPEX_DENSE, SPEX_MPZ, A, option);
    option->print_level = 3;
    SPEX_QR_IPGE( A2, &R2, &Q2);
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


    SPEX_info info;
    SPEX_matrix rhos = NULL,R3=NULL, rhos2 = NULL;
    int64_t *h;
    int64_t j=0, nnz=n*m;
    int64_t i,pQ, pR;
    int sgn;
    SPEX_symbolic_analysis S = NULL;
    SPEX_factorization F = NULL ;

    printf("analysis:\n");
    option->print_level = 3;
    option->order =  SPEX_NO_ORDERING;
    DEMO_OK (SPEX_qr_analyze(&S, A, option));

    printf("facts:\n");
    option->print_level = 3;
    DEMO_OK (SPEX_qr_factorize(&F, A, S, option));
    
    printf("solve:\n");
    DEMO_OK (SPEX_qr_solve(&x, F, b, option));
    printf("orint x sparse:\n");
    SPEX_matrix_check(x, option);
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Declare permuted matrix and S
    /*SPEX_matrix PAQ = NULL, ATA = NULL;
    SPEX_symbolic_analysis S = NULL;
    // Declare local variables for symbolic analysis
    int64_t *post = NULL;
    int64_t *c = NULL, *cInv=NULL;
    int64_t i, nz;
    
    char *mat_name2 = "ExampleMats/smallZerosATA.mat.txt";
    //char *mat_name = "ExampleMats/LF10.mat.txt";
    //char *rhs_name = "ExampleMats/LF10.rhs.txt";
    // Read in A
    FILE *mat_file2 = fopen(mat_name2,"r");
    if( mat_file2 == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }

DEMO_OK(spex_demo_tripread(&ATA, mat_file2, SPEX_FP64, option));
    fclose(mat_file);

    //--------------------------------------------------------------------------
    // Preorder: obtain the row/column ordering of ATA (Default is COLAMD)
    //--------------------------------------------------------------------------
/*
option->order = 0;
    DEMO_OK( spex_cholesky_preorder(&S, ATA, option) );
    S->Q_perm=S->Pinv_perm;
    //--------------------------------------------------------------------------
    // Permute matrix A, that is apply the row/column ordering from the
    // symbolic analysis step to get the permuted matrix PAQ.
    //--------------------------------------------------------------------------

    DEMO_OK( spex_qr_permute_A(&PAQ, A, true, S, option)); //TODO can make false when you can transpose an empty matrix 
    //DEMO_OK( spex_cholesky_permute_A(&PAQ, ATA, true, S)); 
    option->print_level = 3;
    SPEX_matrix_check(PAQ, option);
    
    //--------------------------------------------------------------------------
    // Symbolic Analysis: compute the elimination tree of PAQ
    //--------------------------------------------------------------------------

    // Obtain elimination tree of A
    DEMO_OK( spex_qr_etree(&S->parent, PAQ) );
    

    // Postorder the elimination tree of A
    DEMO_OK( spex_cholesky_post(&post, S->parent, n) );


    // Get the column counts of A
    DEMO_OK( spex_qr_counts(&c, PAQ, S->parent, post) ); //c is S->cp but backwards

    // Set the column pointers of R
    S->cp = (int64_t*) SPEX_malloc( (n+1)*sizeof(int64_t*));
    if (S->cp == NULL)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    cInv = (int64_t*) SPEX_malloc(n* sizeof (int64_t));
    for(i=1;i<=n;i++)
    {
        cInv[i]=c[n-i];//FIXME there has to be a better way of doing this check L vs R
    }
    
    for(i=0;i<=n;i++)
    {
        //cInv[i]=c[n-i];//FIXME there has to be a better way of doing this check L vs R
        printf("%ld cinv %ld c %ld\n",i, cInv[i],c[i]);
    }
    
    
    /*SPEX_factorization F = NULL ;

    printf("analysis:\n");
    option->print_level = 3;
    option->order = SPEX_NO_ORDERING;
    DEMO_OK (SPEX_cholesky_analyze(&S, ATA, option));

    printf("facts:\n");
    option->print_level = 3;
    DEMO_OK (SPEX_cholesky_factorize(&F, ATA, S, option));
    SPEX_matrix_check(F->L, option);
    */
    
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;

}

