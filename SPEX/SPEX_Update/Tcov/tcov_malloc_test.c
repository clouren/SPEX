//------------------------------------------------------------------------------
// SPEX_CHOLMOD/Tcov/tcov_malloc_test.c
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

#include "tcov_malloc_test.h"

int64_t malloc_count = INT64_MAX ;

// Note that only the ANSI C memory manager is used here
// (malloc, calloc, realloc, free)

// wrapper for malloc
void *tcov_malloc
(
    size_t size        // Size to alloc
)
{
    if (--malloc_count < 0)
    {
        /* pretend to fail */
        printf("malloc pretend to fail\n");
        return (NULL) ;
    }
    return (malloc (size)) ;
}

// wrapper for calloc
void *tcov_calloc
(
    size_t n,          // Size of array
    size_t size        // Size to alloc
)
{
    if (--malloc_count < 0)
    {
        /* pretend to fail */
        printf ("calloc pretend to fail\n");
        return (NULL) ;
    }
    // ensure at least one byte is calloc'd
    return (calloc (n, size)) ;
}

// wrapper for realloc
void *tcov_realloc
(
    void *p,           // Pointer to be realloced
    size_t new_size    // Size to alloc
)
{
    if (--malloc_count < 0)
    {
        /* pretend to fail */
        printf("realloc pretend to fail\n");
        return (NULL);
    }
    return (realloc (p, new_size)) ;
}

// wrapper for free
void tcov_free
(
    void *p            // Pointer to be free
)
{
    // This not really needed, but placed here anyway in case the Tcov tests
    // want to do something different that free(p) in the future.
    free (p) ;
}

jmp_buf spex_gmp_environment ;  // for setjmp and longjmp

int spex_gmp_realloc_test
(
    void **p_new,
    void * p_old,
    size_t old_size,
    size_t new_size
)
{
    int spex_gmp_status = setjmp (spex_gmp_environment);
    if (spex_gmp_status != 0)
    {
        return SPEX_OUT_OF_MEMORY;
    }
    *p_new = spex_gmp_reallocate(p_old, old_size, new_size);
    return SPEX_OK;
}
