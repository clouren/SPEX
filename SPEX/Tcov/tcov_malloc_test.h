//------------------------------------------------------------------------------
// SPEX/SPEX/Tcov/tcov_malloc_test.h
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#ifndef TCOV_SPEX_MALLOC_TEST_H
#define TCOV_SPEX_MALLOC_TEST_H

#include "spex_left_lu_internal.h"
#include "spex_update_internal.h"
#include "SPEX_gmp.h"

extern int64_t malloc_count ;

#ifdef GOTCHA
#undef GOTCHA
#endif
#define GOTCHA \
    printf ("%s, line %d, spex_gmp_ntrials = %ld, malloc_count = %ld\n", \
    __FILE__, __LINE__, spex_gmp_ntrials, malloc_count);

#define SPEX_PRINT_INFO(info)                                               \
{                                                                           \
    printf ("file %s line %d: ", __FILE__, __LINE__) ;                  \
    switch(info)                                                            \
    {                                                                       \
        case SPEX_OK:              printf("SPEX_OK\n");            break;   \
        case SPEX_OUT_OF_MEMORY:   printf("OUT OF MEMORY\n");      break;   \
        case SPEX_SINGULAR:        printf("Matrix is SINGULAR\n"); break;   \
        case SPEX_INCORRECT_INPUT: printf("INCORRECT INPUT\n");    break;   \
        case SPEX_INCORRECT:       printf("SPEX_INCORRECT\n");     break;   \
        default:                   printf("unknown!\n");                    \
    }                                                                       \
}

// wrapper for SPEX_initialize*, SPEX_finalize, and all SPEX_*free functions
#define TEST_OK(method)                             \
if (!pretend_to_fail)                               \
{                                                   \
     info = (method) ; assert (info == SPEX_OK) ;   \
}

//wrapper for all other SPEX_* functions
#define TEST_CHECK(method)                       \
if (!pretend_to_fail)                            \
{                                                \
    info = (method) ;                            \
    if (info == SPEX_OUT_OF_MEMORY)              \
    {                                            \
        SPEX_FREE_ALL;                           \
        pretend_to_fail = true ;                 \
    }                                            \
    else if (info != SPEX_OK)                    \
    {                                            \
        SPEX_PRINT_INFO (info) ;                 \
        printf ("test failure at line %d\n", __LINE__) ;         \
        abort ( ) ;                              \
    }                                            \
}

// wrapper for SPEX_* function when expected error would produce
#define TEST_CHECK_FAILURE(method,expected_error)               \
if (!pretend_to_fail)                            \
{                                                \
    info = (method) ;                            \
    if (info == SPEX_OUT_OF_MEMORY)              \
    {                                            \
        SPEX_FREE_ALL;                           \
        pretend_to_fail = true ;                 \
    }                                            \
    else if (info != expected_error)             \
    {                                            \
        printf ("SPEX method was expected to fail, but succeeded!\n") ; \
        printf ("this error was expected:\n") ;  \
        SPEX_PRINT_INFO (expected_error) ;       \
        printf ("but this error was obtained:\n") ;     \
        SPEX_PRINT_INFO (info) ;                 \
        printf ("test failure at line %d (wrong error)\n", __LINE__) ;  \
        abort ( ) ;                              \
    }                                            \
}

// wrapper for malloc
void *tcov_malloc
(
    size_t size        // Size to alloc
) ;

// wrapper for calloc
void *tcov_calloc
(
    size_t n,          // Size of array
    size_t size        // Size to alloc
) ;

// wrapper for realloc
void *tcov_realloc
(
    void *p,           // Pointer to be realloced
    size_t new_size    // Size to alloc
) ;

// wrapper for free
void tcov_free
(
    void *p            // Pointer to be free
) ;

// used to test spex_gmp_reallocate
int spex_gmp_realloc_test
(
    void **p_new,
    void * p_old,
    size_t old_size,
    size_t new_size
);

SPEX_info spex_check_solution
(
    bool *Is_correct,             // if the solution is correct
    const SPEX_matrix *A,         // Input matrix of CSC MPZ
    const SPEX_matrix *x,         // Solution vectors
    const SPEX_matrix *b,         // Right hand side vectors
    const SPEX_options* option    // Command options
);

SPEX_info spex_update_verify
(
    bool *Is_correct,         // if the factorization is correct
    SPEX_factorization *F,    // LU factorization of A
    const SPEX_matrix *A,     // Input matrix of SPEX_DYNAMIC_CSC MPZ
    const SPEX_options *option// command options
);
#endif

