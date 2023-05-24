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
#include "spex_demos.h"

#define FREE_WORKSPACE                          \
    SPEX_matrix_free(&A,  NULL);                \
    SPEX_matrix_free(&A2, NULL);                \
    SPEX_matrix_free(&R,  NULL);                \
    SPEX_matrix_free(&Q,  NULL);                \
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

//m=5;n=5;seed=14;
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
    SPEX_matrix Q = NULL;
    SPEX_matrix R = NULL;
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

    SPEX_generate_random_matrix ( &Ainit, m, n, seed, lower, upper);
    Ainit->nz = m*n;

    // Create A as a copy of Ainit
    // A is a copy of the Ainit matrix. A is a sparse matrix with mpz_t entries
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, Ainit, option);

    // Create A2 as a copy of Ainit
    // A2 is a copy of the Ainit matrix. A is a dense matrix with mpz_t entries
    SPEX_matrix_copy(&A2, SPEX_DENSE, SPEX_MPZ, Ainit, option);

    //--------------------------------------------------------------------------
    // Factorize
    //--------------------------------------------------------------------------


    SPEX_QR_IPGE( A2, &R2, &Q2);

    // playing with sparse
    //will need to change DEMO_OK back to DEMO_OK
    SPEX_info info;
    SPEX_matrix rhos = NULL,R3=NULL, rhos2 = NULL;
    int64_t *h;
    int64_t j=0, nnz=n*m;
    int64_t i,pQ, pR;
    int sgn;

    h = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    DEMO_OK (SPEX_matrix_allocate(&(rhos2), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));

    // Allocate R. We are performing the Thin REF QR factorization so
    // R is n*n
    //DEMO_OK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, n*n, false, true, NULL));
    //SPEX_matrix_copy(&R, SPEX_CSC, SPEX_MPZ, R2, option);

    //Q=A (will need to change as soon as A is actually sparse)
    DEMO_OK(SPEX_matrix_copy(&Q, SPEX_CSC, SPEX_MPZ, A, NULL));
    //will also need to change once A is actually sparse (and we're actually only doing one col)
    for (i = 0; i < nnz; i++)
    {
        h[i] = -1;
    }

    for (i = 0; i < n; i++)
    {
        SPEX_MPZ_SET(rhos2->x.mpz[i], SPEX_2D(R2, i, i, mpz)); //rho^i=R2(i,i)
    }

    option->print_level = 3;
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, Ainit, option);
    /*SPEX_matrix_check(A, option); //checa que onda, a lo mejor es lo de option
    printf("finishA\n");
    SPEX_matrix_check(R2, option);
    SPEX_matrix_check(Q2, option);*/

    
    //R
    SPEX_symbolic_analysis S = NULL;
    /*SPEX_matrix PAQ = NULL;
    int64_t *post = NULL;
    int64_t *c = NULL;

    DEMO_OK( spex_qr_preorder(&S, A, option) );

    DEMO_OK( spex_qr_permute_A(&PAQ, A, false, S,option) );

    // Obtain elimination tree of A
    DEMO_OK( spex_qr_etree(&S->parent, A) );

    // Postorder the elimination tree of A
    DEMO_OK( spex_cholesky_post(&post, S->parent, n) );

    // Get the column counts of A
    DEMO_OK( spex_qr_counts(&c, A, S->parent, post) );*/
    /*DEMO_OK (SPEX_qr_analyze(&S, A, option));
   

    SPEX_CHECK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, S->unz,
                                    false, false, option));
    printf("here\n");
    // Set the column pointers of L
    for (int64_t k = 0; k < n; k++)
    {   
        R->p[k] = S->cp[k];
        printf("k %ld, p[k] %ld\n", k, R->p[k]);
    }
*/ /*printf("sparse\n");
    SPEX_matrix_check(R, option);
    SPEX_matrix_check(Q, option);
    printf("dense\n");
    SPEX_matrix_check(R2, option);
    SPEX_matrix_check(Q2, option);*/
    SPEX_factorization F = NULL ;

    printf("before analysis\n");
    DEMO_OK (SPEX_qr_analyze(&S, A, option));
   printf("after analysis\n");
    //DEMO_OK (SPEX_qr_factorize(&R, &Q, &rhos, A, S, option));
    DEMO_OK (SPEX_qr_factorize(&F, A, S, option));
    printf("after fact\n");
   
    SPEX_generate_random_matrix ( &b2, m, 1, seed, lower, upper);
    b2->nz = m;
    // Make a copy of b
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b2, option);

    DEMO_OK (SPEX_qr_solve(&x, F, b, option));
    printf("orint x sparse:\n");
    SPEX_matrix_check(x, option);

    SPEX_Qtb(Q2, b, &b_new);
    //spex_matrix_mul(b_new,R->x.mpz[R->nz]);

    SPEX_QR_backsolve(R2, b_new, &x2);
    printf("orint x dense:\n");
    SPEX_matrix_check(x2, option);
    /*SPEX_mpq_set_num(x->scale, SPEX_2D(R, n-1, n-1, mpz));


    // Create a double version of x
    SPEX_matrix_copy(&x_doub, SPEX_DENSE, SPEX_FP64, x, option);
    SPEX_matrix_check(b_new, option);*/

   /* DEMO_OK(SPEX_matrix_allocate(&R, SPEX_CSC, SPEX_MPZ, n, n, S->unz,
                                    false, false, option));
    // Set the column pointers of R
    for (int64_t k = 0; k < n; k++)
    {   
        R->p[k] = S->cp[k+1];
        printf("k %ld, p[k] %ld\n", k, R->p[k]);
    }
    R->p[n] = S->unz;

    R->i[0]=0;
    R->i[1]=0;
    R->i[2]=1;
    R->i[3]=0;

    R->i[4]=1;
    R->i[5]=2;
    R->i[6]=0;
    R->i[7]=1;

    R->i[8]=2;
    R->i[9]=3;



*/
  /* spex_qr_pre_factor(&R,A,S);
   printf("here %ld\n",R->nz);
   for (i = 0; i < R->n; i++)
    {
        printf("p: %ld %ld\n",R->p[i],i);
    }
    
    for (i = 0; i < R->nzmax; i++)
    {
        SPEX_MPZ_INIT(R->x.mpz[i]);
        SPEX_MPZ_SET_UI(R->x.mpz[i],i);
    }

    SPEX_matrix_check(R, option);*/
    /*SPEX_matrix_check(Q, option);*/
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;

}

