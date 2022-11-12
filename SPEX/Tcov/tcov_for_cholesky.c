//------------------------------------------------------------------------------
// SPEX/Tcov/tcov_for_cholesky.c: test coverage for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* This program will exactly solve the sparse linear system Ax = b by performing
 * the SPEX Cholesky factorization.
 */

#include "tcov_malloc_test.h"
#include "chol_demos.h"

#include <float.h>
#include <assert.h>

//------------------------------------------------------------------------------
// OK: call a method and assert that it succeeds
//------------------------------------------------------------------------------

// The method must return a SPEX_info value.

#define OK(method)                              \
{                                               \
    info = (method) ;                           \
    if (info != SPEX_OK)                        \
    {                                           \
        printf ("SPEX Cholesky test failure, line %d, info %d\n", \
            __LINE__, info) ;                   \
        fprintf (stderr, "SPEX Cholesky test failure, line %d, info %d\n", \
            __LINE__, info) ;                   \
        abort ( ) ;                             \
    }                                           \
}

//------------------------------------------------------------------------------
// BRUTAL: test a method with debug malloc, until it succeeds
//------------------------------------------------------------------------------

// The method must return a bool (true if successful, false if failure).

#define NTRIAL_MAX 10000

#define BRUTAL(method)                                                      \
{                                                                           \
    bool done = false ;                                                     \
    int64_t trial = 0 ;                                                     \
    for (trial = 0 ; !done && trial <= NTRIAL_MAX ; trial++)                \
    {                                                                       \
        /* printf (":") ; fflush (stdout) ; */                              \
        malloc_count = trial ;                                              \
        done = (method) ;                                                   \
    }                                                                       \
    OK (done ? SPEX_OK : SPEX_PANIC) ;                                      \
    if (done) printf ("\nCholesky trials %ld: all tests passed\n", trial) ; \
}

//------------------------------------------------------------------------------
// read_test_matrix: read in a matrix from a file
//------------------------------------------------------------------------------

void read_test_matrix (SPEX_matrix *A_handle, char *filename) ;

void read_test_matrix (SPEX_matrix *A_handle, char *filename)
{
    SPEX_info info ;
    FILE *f = fopen (filename, "r") ;
    OK (f == NULL ? SPEX_PANIC : SPEX_OK) ;
    OK (SPEX_tripread_double (A_handle, f, NULL)) ;
    fclose (f) ;
}

//------------------------------------------------------------------------------
// create_test_rhs: create a right-hand-side vector
//------------------------------------------------------------------------------

void create_test_rhs (SPEX_matrix *b_handle, int64_t n) ;

void create_test_rhs (SPEX_matrix *b_handle, int64_t n)
{
    SPEX_info info ;
    OK (SPEX_matrix_allocate (b_handle, SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, NULL)) ;
    SPEX_matrix b = *(b_handle) ;
    // b(0)=0
    OK (SPEX_mpz_set_ui (b->x.mpz [0], 0)) ;
    for (int64_t k = 1 ; k < n ; k++)
    {
        // b(k) = 1
        OK (SPEX_mpz_set_ui (b->x.mpz [k], 1)) ;
    }
}

//------------------------------------------------------------------------------
// spex_test_chol_backslash: test SPEX_cholesky_backslash
//------------------------------------------------------------------------------

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                           \
{                                               \
    OK (SPEX_matrix_free (&x, option)) ;        \
}

bool spex_test_chol_backslash (SPEX_matrix A, SPEX_matrix b,
    SPEX_options option) ;

bool spex_test_chol_backslash (SPEX_matrix A, SPEX_matrix b,
    SPEX_options option)
{
    SPEX_info info ;
    bool pretend_to_fail = false ;
    SPEX_matrix x = NULL ;
    // solve Ax=b
    TEST_CHECK (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option)) ;
    // disable memory testing when checking the solution
    int64_t save = malloc_count ; malloc_count = INT64_MAX ;
    TEST_CHECK (SPEX_check_solution (A, x, b, option)) ;
    // re-enable memory testing
    malloc_count = save ;
    SPEX_FREE_ALL ;
    return (!pretend_to_fail) ;
}

//------------------------------------------------------------------------------
// spex_test_chol_afs: test SPEX_cholesky_[analyze,factorize,solve]
//------------------------------------------------------------------------------

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                                   \
{                                                       \
    OK (SPEX_symbolic_analysis_free (&S, option)) ;     \
    OK (SPEX_factorization_free (&F, option)) ;         \
    OK (SPEX_matrix_free (&x, option)) ;                \
}

bool spex_test_chol_afs (SPEX_matrix A, SPEX_matrix b, SPEX_options option) ;

bool spex_test_chol_afs (SPEX_matrix A, SPEX_matrix b, SPEX_options option)
{
    SPEX_info info ;
    bool pretend_to_fail = false ;
    SPEX_symbolic_analysis *S = NULL ;
    SPEX_factorization F = NULL ;
    SPEX_matrix x = NULL ;
    // solve Ax=b
    TEST_CHECK (SPEX_cholesky_analyze (&S, A, option)) ;
    TEST_CHECK (SPEX_cholesky_factorize (&F, A, S, option)) ;
    TEST_CHECK (SPEX_cholesky_solve (&x, F, b, option)) ;
    // disable memory testing when checking the solution
    int64_t save = malloc_count ; malloc_count = INT64_MAX ;
    TEST_CHECK (SPEX_check_solution (A, x, b, option)) ;
    // re-enable memory testing
    malloc_count = save ;
    SPEX_FREE_ALL ;
    return (!pretend_to_fail) ;
}

//------------------------------------------------------------------------------
// tcov_for_cholesky: main program
//------------------------------------------------------------------------------

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                                   \
{                                                       \
    OK (SPEX_symbolic_analysis_free (&S, option)) ;     \
    OK (SPEX_factorization_free (&F, option)) ;         \
    OK (SPEX_matrix_free (&x, option)) ;                \
    OK (SPEX_matrix_free (&A, option)) ;                \
    OK (SPEX_matrix_free (&b, option)) ;                \
    SPEX_FREE (option) ;                                \
}

int main (int argc, char *argv [])
{

    //--------------------------------------------------------------------------
    // start SPEX
    //--------------------------------------------------------------------------

    SPEX_matrix A = NULL, b = NULL, x = NULL ;
    SPEX_symbolic_analysis *S = NULL ;
    SPEX_factorization F = NULL ;
    SPEX_options option = NULL ;
    bool pretend_to_fail = false ;

    if (argc < 2)
    {
        printf ("usage: tcov_for_cholesky matrixfilename\n") ;
        abort ( ) ;
    }

    SPEX_info info ;
    OK (SPEX_initialize_expert (tcov_malloc, tcov_calloc, tcov_realloc,
        tcov_free)) ;

    // disable malloc testing for the first part of the test
    spex_gmp_ntrials = INT64_MAX ;
    malloc_count = INT64_MAX ;

    OK (SPEX_create_default_options (&option)) ;

    //--------------------------------------------------------------------------
    // test a few small invalid matrices
    //--------------------------------------------------------------------------

    // unsymmetric matrix (row counts != col counts)
    printf ("Cholesky: error handling for unsymmetric matrix (1)\n") ;
    read_test_matrix (&A, "../ExampleMats/test1.mat.txt") ;
    create_test_rhs (&b, A->n) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_UNSYMMETRIC) ;
    OK (SPEX_matrix_free (&A, option)) ;
    OK (SPEX_matrix_free (&b, option)) ;

    // unsymmetric matrix (unsymmetric pattern)
    printf ("Cholesky: error handling for unsymmetric matrix (2)\n") ;
    read_test_matrix (&A, "../ExampleMats/test2.mat.txt") ;
    create_test_rhs (&b, A->n) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_UNSYMMETRIC) ;
    OK (SPEX_matrix_free (&A, option)) ;
    OK (SPEX_matrix_free (&b, option)) ;

    // unsymmetric matrix (unsymmetric values)
    printf ("Cholesky: error handling for unsymmetric matrix (3)\n") ;
    read_test_matrix (&A, "../ExampleMats/test3.mat.txt") ;
    create_test_rhs (&b, A->n) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_UNSYMMETRIC) ;
    OK (SPEX_matrix_free (&A, option)) ;
    OK (SPEX_matrix_free (&b, option)) ;

    // symmetric indefinite matrix
    printf ("Cholesky: error handling for symmetric indefinite matrix (4)\n") ;
    read_test_matrix (&A, "../ExampleMats/test4.mat.txt") ;
    create_test_rhs (&b, A->n) ;
    option->algo = SPEX_CHOL_UP ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_NOTSPD) ;
    option->algo = SPEX_CHOL_LEFT ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_NOTSPD) ;
    option->algo = SPEX_CHOL_UP ;
    OK (SPEX_matrix_free (&A, option)) ;
    OK (SPEX_matrix_free (&b, option)) ;

    //--------------------------------------------------------------------------
    // load the test matrix and create the right-hand-side
    //--------------------------------------------------------------------------

    read_test_matrix (&A, argv [1]) ;
    int64_t n = A->n ;
    int64_t m = A->m ;
    int64_t anz = -1 ;
    OK (SPEX_matrix_nnz (&anz, A, option)) ;
    printf ("\nInput matrix: %ld-by-%ld with %ld entries\n", n, m, anz) ;
    OK ((n != m) ? SPEX_PANIC : SPEX_OK) ;
    create_test_rhs (&b, A->n) ;

    //--------------------------------------------------------------------------
    // test SPEX_transpose
    //--------------------------------------------------------------------------

    printf("\n Test SPEX_transpose \n");
    SPEX_matrix A_mpq = NULL, A_mpfr = NULL, A_int = NULL, A_fp = NULL;
    SPEX_matrix T_mpq = NULL, T_mpfr = NULL, T_int = NULL, T_fp = NULL;
    // T = A'
    OK ( SPEX_matrix_copy(&A_mpq, SPEX_CSC, SPEX_MPQ, A, option));
    OK ( SPEX_transpose(&T_mpq, A_mpq, option) );

    OK ( SPEX_matrix_copy(&A_mpfr, SPEX_CSC, SPEX_MPFR, A, option));
    OK ( SPEX_transpose(&T_mpfr, A_mpfr, option) );

    OK ( SPEX_matrix_copy(&A_int, SPEX_CSC, SPEX_INT64, A, option));
    OK ( SPEX_transpose(&T_int, A_int, option) );

    OK ( SPEX_matrix_copy(&A_fp, SPEX_CSC, SPEX_FP64, A, option));
    OK ( SPEX_transpose(&T_fp, A_fp, option));

    SPEX_matrix_free(&A_mpq,option);
    SPEX_matrix_free(&A_mpfr,option);
    SPEX_matrix_free(&A_int,option);
    SPEX_matrix_free(&A_fp,option);
    SPEX_matrix_free(&T_mpq,option);
    SPEX_matrix_free(&T_mpfr,option);
    SPEX_matrix_free(&T_int,option);
    SPEX_matrix_free(&T_fp,option);
    

    //--------------------------------------------------------------------------
    // error handling
    //--------------------------------------------------------------------------

    // inputs cannot be NULL
    TEST_CHECK_FAILURE (SPEX_matrix_nnz (NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_matrix_nnz (NULL, A, NULL),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_matrix_nnz (&anz, NULL, NULL),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_analyze (NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (NULL, SPEX_MPQ, NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_factorize (NULL, NULL, NULL, NULL),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_determine_symmetry (NULL, NULL ),
        SPEX_INCORRECT_INPUT) ;

    // type cannot be int64
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_INT64, A, b, option),
        SPEX_INCORRECT_INPUT) ;

    // mangle the matrix: invalid dimensions
    A->n = 0 ;
    A->m = 0 ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_INCORRECT_INPUT) ;
    A->n = n ;
    A->m = n ;

    // mangle the matrix: invalid type
    A->type = SPEX_INT64 ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_INCORRECT_INPUT) ;
    TEST_CHECK_FAILURE (SPEX_determine_symmetry (A, option),
        SPEX_INCORRECT_INPUT) ;
    A->type = SPEX_MPZ ;

    // valid analysis, but break the factorization
    OK (SPEX_cholesky_analyze (&S, A, option)) ;
    A->type = SPEX_INT64 ;
    TEST_CHECK_FAILURE (SPEX_cholesky_factorize (&F, A, S, option),
        SPEX_INCORRECT_INPUT) ;
    A->type = SPEX_MPZ ;
    OK (SPEX_symbolic_analysis_free (&S, option)) ;

    // valid analysis and factorization, but break the solve
    OK (SPEX_cholesky_analyze (&S, A, option)) ;
    OK (SPEX_cholesky_factorize (&F, A, S, option)) ;
    b->type = SPEX_INT64 ;
    TEST_CHECK_FAILURE (SPEX_cholesky_solve (&x, F, b, option),
        SPEX_INCORRECT_INPUT) ;
    b->type = SPEX_MPZ ;
    OK (SPEX_symbolic_analysis_free (&S, option)) ;
    OK (SPEX_factorization_free (&F, option)) ;

    // invalid algorithm
    option->algo = SPEX_QR_GRAM ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option),
        SPEX_INCORRECT_ALGORITHM) ;
    option->algo = SPEX_CHOL_UP ;

    //--------------------------------------------------------------------------
    // solve Ax=b with SPEX_cholesky_backslash and check the solution
    //--------------------------------------------------------------------------

    option->order = SPEX_AMD ;
    option->algo = SPEX_CHOL_UP ;
    option->print_level = 3 ;
    printf ("Cholesky backslash, up-looking, no malloc testing:\n") ;
    bool ok = spex_test_chol_backslash (A, b, option) ;
    OK (ok ? SPEX_OK : SPEX_PANIC) ;
    option->print_level = 0 ;

    printf ("Cholesky backslash, up-looking, no malloc testing, colamd:\n") ;
    option->order = SPEX_COLAMD ;
    option->print_level = 3 ;
    ok = spex_test_chol_backslash (A, b, option) ;
    OK (ok ? SPEX_OK : SPEX_PANIC) ;
    option->order = SPEX_AMD ;
    option->print_level = 0 ;

    printf ("Cholesky backslash, no malloc testing, natural ordering:\n") ;
    option->order = SPEX_NO_ORDERING ;
    ok = spex_test_chol_backslash (A, b, option) ;
    OK (ok ? SPEX_OK : SPEX_PANIC) ;

    printf ("Cholesky backslash, no malloc testing, return x as MPFR:\n") ;
    TEST_CHECK (SPEX_cholesky_backslash (&x, SPEX_MPFR, A, b, option)) ;
    //NOTE: mpfr solution can't be checked because mpfr->mpz isn't guaranteed
    //      to be exact
    OK (SPEX_matrix_free (&x, option)) ;

    printf ("Cholesky backslash, up-looking with malloc testing, colamd:\n") ;
    option->order = SPEX_COLAMD ;
    BRUTAL (spex_test_chol_backslash (A, b, option)) ;

    printf ("Cholesky backslash, up-looking with malloc testing, amd:\n") ;
    option->order = SPEX_AMD ;
    BRUTAL (spex_test_chol_backslash (A, b, option)) ;
    
    printf ("Cholesky backslash, up-looking with malloc testing, no ordering:\n") ;
    option->order = SPEX_NO_ORDERING ;
    BRUTAL (spex_test_chol_backslash (A, b, option)) ;

    printf ("Cholesky backslash, left-looking with malloc testing:\n") ;
    option->algo = SPEX_CHOL_LEFT ;
    BRUTAL (spex_test_chol_backslash (A, b, option)) ;

    //--------------------------------------------------------------------------
    // solve Ax=b with SPEX_cholesky_[analyze,factorize,solve] and check solution
    //--------------------------------------------------------------------------

    option->algo = SPEX_CHOL_UP ;

    printf ("Cholesky analyze/factorize/solve, no malloc testing:\n") ;
    spex_gmp_ntrials = INT64_MAX ;
    malloc_count = INT64_MAX ;
    ok = spex_test_chol_afs (A, b, option) ;
    OK (ok ? SPEX_OK : SPEX_PANIC) ;

    printf ("Cholesky analyze/factorize/solve, with malloc testing:\n") ;
    // also check a different RHS, with b(0) = 0
    OK (SPEX_mpz_set_ui (b->x.mpz [0], 0)) ;
    BRUTAL (spex_test_chol_afs (A, b, option)) ;

    //--------------------------------------------------------------------------
    // free the test problem and finalize SPEX
    //--------------------------------------------------------------------------

    SPEX_FREE_ALL ;
    OK (SPEX_finalize ( )) ;

    //--------------------------------------------------------------------------
    // error handling
    //--------------------------------------------------------------------------

    // SPEX not initialized
    TEST_CHECK_FAILURE (SPEX_cholesky_analyze (NULL, NULL, NULL), SPEX_PANIC) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_solve (NULL, NULL, NULL, NULL), SPEX_PANIC) ;
    TEST_CHECK_FAILURE (SPEX_cholesky_backslash (NULL, SPEX_MPQ, NULL, NULL, NULL),
        SPEX_PANIC) ;

    printf ("\nSPEX_cholesky: all tests passed\n") ;
    fprintf (stderr, "%s: all tests passed\n\n", __FILE__) ;
}

