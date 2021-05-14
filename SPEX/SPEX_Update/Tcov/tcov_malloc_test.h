//------------------------------------------------------------------------------
// SPEX_CHOLMOD/Tcov/tcov_malloc_test.h
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

#ifndef TCOV_SPEX_MALLOC_TEST_H
#define TCOV_SPEX_MALLOC_TEST_H

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
    printf ("file %s line %d: ", __FILE__, __LINE__) ;                      \
    switch(info)                                                            \
    {                                                                       \
        case SPEX_OK:              printf("SPEX_OK\n");            break;   \
        case SPEX_OUT_OF_MEMORY:   printf("OUT OF MEMORY\n");      break;   \
        case SPEX_SINGULAR:        printf("Matrix is SINGULAR\n"); break;   \
        case SPEX_INCORRECT_INPUT: printf("INCORRECT INPUT\n");    break;   \
        default:                   printf("unknown!\n");                    \
    }                                                                       \
}

#ifdef SPEX_CHECK
#undef SPEX_CHECK
#endif

#define SPEX_CHECK(method)          \
{                                   \
    info = (method) ;               \
    if (info != SPEX_OK)            \
    {                               \
        SPEX_PRINT_INFO (info)      \
        SPEX_FREE_ALL ;             \
        return (info) ;             \
    }                               \
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
#endif

