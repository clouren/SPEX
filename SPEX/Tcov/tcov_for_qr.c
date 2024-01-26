// ----------------------------------------------------------------------------
// SPEX/Tcov/tcov_for_qr.c: test coverage for SPEX_QR
// ----------------------------------------------------------------------------

// SPEX: (c) 2019-2023, Chris Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//-----------------------------------------------------------------------------

/* This program will exactly solve the sparse linear system Ax = b by performing
 * the SPEX QR factorization.
 */

#include "tcov_utilities.h"
#include "spex_demos.h"

// test wrapper for SPEX_* function when expected error would produce
#define ERR(method,expected_error)                                      \
{                                                                       \
    SPEX_info info5 = (method) ;                                        \
    if (info5 != expected_error)                                        \
    {                                                                   \
        printf ("SPEX method was expected to fail, but succeeded!\n") ; \
        printf ("this error was expected:\n") ;                         \
        SPEX_PRINT_INFO (expected_error) ;                              \
        printf ("but this error was obtained:\n") ;                     \
        TEST_ABORT (info5) ;                                            \
    }                                                                   \
}

//------------------------------------------------------------------------------
// BRUTAL: test a method with debug malloc, until it succeeds
//------------------------------------------------------------------------------

// The method must return a bool (true if successful, false if failure).

#define NTRIAL_MAX 1//10000

#define BRUTAL(method)                                                      \
{                                                                           \
    int64_t trial = 0 ;                                                     \
    SPEX_info info2 = SPEX_OK ;                                             \
    for (trial = 0 ; trial <= NTRIAL_MAX ; trial++)                         \
    {                                                                       \
        malloc_count = trial ;                                              \
        info2 = (method) ;                                                  \
        if (info2 != SPEX_OUT_OF_MEMORY) break ;                            \
    }                                                                       \
    if (info2 != SPEX_OK) TEST_ABORT (info2) ;                              \
    malloc_count = UINT64_MAX ;                                             \
    printf ("\nBrutal QR trials %ld: tests passed\n", trial);         \
}

//------------------------------------------------------------------------------
// read_test_matrix: read in a matrix from a file
//------------------------------------------------------------------------------

void read_test_matrix (SPEX_matrix *A_handle, char *filename);

void read_test_matrix (SPEX_matrix *A_handle, char *filename)
{
    FILE *f = fopen (filename, "r");
    OK (f == NULL ? SPEX_PANIC : SPEX_OK);
    OK (spex_demo_tripread (A_handle, f, SPEX_FP64, NULL));
    fclose (f);
}

//------------------------------------------------------------------------------
// create_test_rhs: create a right-hand-side vector
//------------------------------------------------------------------------------

void create_test_rhs (SPEX_matrix *b_handle, int64_t n);

void create_test_rhs (SPEX_matrix *b_handle, int64_t n)
{
    OK (SPEX_matrix_allocate (b_handle, SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, NULL));
    SPEX_matrix b = *(b_handle);
    // b(0)=0
    OK (SPEX_mpz_set_ui (b->x.mpz [0], 0));
    for (int64_t k = 1 ; k < n ; k++)
    {
        // b(k) = 1
        OK (SPEX_mpz_set_ui (b->x.mpz [k], 1));
    }
}

//------------------------------------------------------------------------------
// spex_test_chol_backslash: test SPEX_qr_backslash
//------------------------------------------------------------------------------

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                           \
{                                               \
    OK (SPEX_matrix_free (&x, option));         \
}

SPEX_info spex_test_qr_backslash (SPEX_matrix A, SPEX_matrix b,
    SPEX_options option);

SPEX_info spex_test_qr_backslash (SPEX_matrix A, SPEX_matrix b,
    SPEX_options option)
{
    SPEX_matrix x = NULL ;
    // solve Ax=b
    OK2 (SPEX_qr_backslash (&x, SPEX_MPQ, A, b, option));
    // disable memory testing when checking the solution
    int64_t save = malloc_count ; malloc_count = INT64_MAX ;
    OK (spex_demo_check_solution (A, x, b, option));
    // re-enable memory testing
    malloc_count = save ;
    SPEX_FREE_ALL;
    return (SPEX_OK) ;
}


//------------------------------------------------------------------------------
// spex_test_chol_afs: test SPEX_qr_[analyze,factorize,solve]
//------------------------------------------------------------------------------

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                                   \
{                                                       \
    OK (SPEX_symbolic_analysis_free (&S, option));      \
    OK (SPEX_factorization_free (&F, option));          \
    OK (SPEX_matrix_free (&x, option));                 \
}

SPEX_info spex_test_qr_afs
(
    SPEX_matrix A,
    SPEX_matrix b,
    SPEX_options option
) ;

SPEX_info spex_test_qr_afs (SPEX_matrix A, SPEX_matrix b, SPEX_options option)
{
    SPEX_symbolic_analysis S = NULL ;
    SPEX_factorization F = NULL ;
    SPEX_matrix x = NULL ;
    // solve Ax=b
    OK2 (SPEX_qr_analyze (&S, A, option));
    OK2 (SPEX_qr_factorize (&F, A, S, option));
    OK2 (SPEX_qr_solve (&x, F, b, option));
    // disable memory testing when checking the solution
    int64_t save = malloc_count ; malloc_count = INT64_MAX ;
    OK (spex_demo_check_solution (A, x, b, option));
    // re-enable memory testing
    malloc_count = save ;
    SPEX_FREE_ALL;
    return (SPEX_OK);
}

//------------------------------------------------------------------------------
// tcov_for_qr: main program
//------------------------------------------------------------------------------

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                                   \
{                                                       \
    OK (SPEX_symbolic_analysis_free (&S, option));      \
    OK (SPEX_factorization_free (&F, option));          \
    OK (SPEX_matrix_free (&x, option));                 \
    OK (SPEX_matrix_free (&A, option));                 \
    OK (SPEX_matrix_free (&b, option));                 \
    SPEX_FREE (option);                                 \
}

int main (int argc, char *argv [])
{

    //--------------------------------------------------------------------------
    // start SPEX
    //--------------------------------------------------------------------------

    SPEX_matrix A = NULL, b = NULL, x = NULL ;
    SPEX_symbolic_analysis S = NULL ;
    SPEX_factorization F = NULL, F2 = NULL ;
    SPEX_options option = NULL ;

    if (argc < 2)
    {
        printf ("usage: tcov_for_qr matrixfilename\n");
        TEST_ABORT (SPEX_INCORRECT_INPUT);
    }

    SPEX_info info ;
    OK (SPEX_initialize_expert (tcov_malloc, tcov_calloc, tcov_realloc,
        tcov_free));

    // disable malloc testing for the first part of the test
    spex_set_gmp_ntrials (INT64_MAX) ;
    malloc_count = INT64_MAX ;

    OK (SPEX_create_default_options (&option));

    //--------------------------------------------------------------------------
    // test a few small invalid matrices
    //--------------------------------------------------------------------------

    // wrong shape matrix (row counts <= col counts)
   /* printf ("QR: error handling for m<n matrix\n");
    read_test_matrix (&A, "../ExampleMats/test1.mat.txt"); //TODO
    create_test_rhs (&b, A->n);
    ERR (SPEX_qr_backslash (&x, SPEX_MPQ, A, b, option), SPEX_NOTSPD);
    OK (SPEX_matrix_free (&A, option));
    OK (SPEX_matrix_free (&b, option));*/

    
    //--------------------------------------------------------------------------
    // load the test matrix and create the right-hand-side
    //--------------------------------------------------------------------------

    read_test_matrix (&A, argv [1]);
    int64_t n = A->n ;
    int64_t m = A->m ;
    int64_t anz = -1 ;
    OK (SPEX_matrix_nnz (&anz, A, option));
    printf ("\nInput matrix: %ld-by-%ld with %ld entries\n", n, m, anz);
    OK ((n != m) ? SPEX_PANIC : SPEX_OK);
    create_test_rhs (&b, A->n);

    //--------------------------------------------------------------------------
    // error handling
    //--------------------------------------------------------------------------

    // inputs cannot be NULL
    ERR (SPEX_matrix_nnz (NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT);
    ERR (SPEX_matrix_nnz (NULL, A, NULL),
        SPEX_INCORRECT_INPUT);
    ERR (SPEX_matrix_nnz (&anz, NULL, NULL),
        SPEX_INCORRECT_INPUT);
    ERR (SPEX_qr_analyze (NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT);
    ERR (SPEX_qr_backslash (NULL, SPEX_MPQ, NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT);
    ERR (SPEX_qr_factorize (NULL, NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT);

    // type cannot be int64
    ERR (SPEX_qr_backslash (&x, SPEX_INT64, A, b, option),
        SPEX_INCORRECT_INPUT);

    // mangle the matrix: invalid dimensions
    A->n = 0 ;
    A->m = 0 ;
    ERR (SPEX_qr_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_INCORRECT_INPUT);
    A->n = n ;
    A->m = m ;

    // mangle the matrix: invalid type
    A->type = SPEX_INT64 ;
    ERR (SPEX_qr_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_INCORRECT_INPUT);

    // valid analysis, but break the factorization
    OK (SPEX_qr_analyze (&S, A, option));
    A->type = SPEX_INT64 ;
    ERR (SPEX_qr_factorize (&F, A, S, option),
        SPEX_INCORRECT_INPUT);
    A->type = SPEX_MPZ ;
    OK (SPEX_symbolic_analysis_free (&S, option));

    // valid analysis and factorization, but break the solve
    OK (SPEX_qr_analyze (&S, A, option));
    OK (SPEX_qr_factorize (&F, A, S, option));
    b->type = SPEX_INT64 ;
    ERR (SPEX_qr_solve (&x, F, b, option),
        SPEX_INCORRECT_INPUT);
    b->type = SPEX_MPZ ;
    OK (SPEX_symbolic_analysis_free (&S, option));
    OK (SPEX_factorization_free (&F, option));

    // invalid algorithm
    option->algo = 99 ;
    ERR (SPEX_qr_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_INCORRECT_ALGORITHM);
    option->algo = SPEX_CHOL_UP ;

    //--------------------------------------------------------------------------
    // solve Ax=b with SPEX_qr_backslash and check the solution
    //--------------------------------------------------------------------------

    option->order = SPEX_COLAMD ;
    option->print_level = 3 ;
    printf ("QR backslash, no malloc testing:\n");
    OK (spex_test_qr_backslash (A, b, option));
    option->print_level = 0 ;

    printf ("QR backslash, no malloc testing, amd:\n");
    option->order = SPEX_AMD ;
    option->print_level = 3 ;
    OK (spex_test_qr_backslash (A, b, option));
    option->order = SPEX_AMD ;
    option->print_level = 0 ;

    printf ("QR backslash, no malloc testing, natural ordering:\n");
    option->order = SPEX_NO_ORDERING ;
    OK (spex_test_qr_backslash (A, b, option));

    printf ("QR backslash, no malloc testing, return x as MPFR:\n");
    OK (SPEX_qr_backslash (&x, SPEX_MPFR, A, b, option));
    //NOTE: mpfr solution can't be checked because mpfr->mpz isn't guaranteed
    //      to be exact
    OK (SPEX_matrix_free (&x, option));

    //--------------------------------------------------------------------------
    // solve Ax=b with SPEX_qr_[analyze,factorize,solve]; check solution
    //--------------------------------------------------------------------------

    option->algo = SPEX_CHOL_UP ;

    printf ("QR analyze/factorize/solve, no malloc testing:\n");
    spex_set_gmp_ntrials (INT64_MAX) ;
    malloc_count = INT64_MAX ;
    OK (spex_test_qr_afs (A, b, option));

    printf ("QR analyze/factorize/solve, with malloc testing:\n");
    // also check a different RHS, with b(0) = 0
    OK (SPEX_mpz_set_ui (b->x.mpz [0], 0));
    BRUTAL (spex_test_qr_afs (A, b, option));

    //--------------------------------------------------------------------------
    // error handling
    //--------------------------------------------------------------------------

    // SPEX not initialized
    spex_set_initialized (false);
    ERR (SPEX_qr_factorize (&F2, A, S, option), SPEX_PANIC);
    ERR (SPEX_qr_analyze (NULL, NULL, NULL), SPEX_PANIC);
    ERR (SPEX_qr_solve (NULL, NULL, NULL, NULL), SPEX_PANIC);
    ERR (SPEX_qr_backslash (NULL, SPEX_MPQ, NULL, NULL, NULL),
        SPEX_PANIC);
    spex_set_initialized (true);
    SPEX_FREE_ALL;

    printf ("%s: all tests passed\n\n", __FILE__);
    fprintf (stderr, "%s: all tests passed\n\n", __FILE__);
    return 0;
}

