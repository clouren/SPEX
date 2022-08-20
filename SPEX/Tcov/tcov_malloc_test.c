//------------------------------------------------------------------------------
// SPEX/SPEX/Tcov/tcov_malloc_test.c
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

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
//      printf("malloc pretend to fail\n");
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
//      printf ("calloc pretend to fail\n");
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
//      printf("realloc pretend to fail\n");
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

extern jmp_buf spex_gmp_environment ;  // for setjmp and longjmp

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

//------------------------------------------------------------------------------
// spex_check_solution: check solution to Ax=b
//------------------------------------------------------------------------------

/* Purpose: Given a solution vector x, check the solution of the linear system
 * Ax = b. This is done by computing a rational-arthmetic A*x == b. This
 * function is provided here only used for debugging purposes, as the routines
 * within SPEX are gauranteed to be exact.
 */

#undef SPEX_FREE_ALL
#define SPEX_FREE_ALL                       \
    SPEX_MPQ_CLEAR(temp);                   \
    SPEX_MPQ_CLEAR(scale);                  \
    SPEX_matrix_free(&b2, NULL);

SPEX_info spex_check_solution
(
    bool *Is_correct,             // if the solution is correct
    const SPEX_matrix *A,         // Input matrix of CSC MPZ
    const SPEX_matrix *x,         // Solution vectors
    const SPEX_matrix *b,         // Right hand side vectors
    const SPEX_options* option    // Command options
)
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    //--------------------------------------------------------------------------
    // check inputs. Input are also checked by the two callers
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_REQUIRE (A, SPEX_CSC,   SPEX_MPZ) ;
    SPEX_REQUIRE (x, SPEX_DENSE, SPEX_MPQ) ;
    SPEX_REQUIRE (b, SPEX_DENSE, SPEX_MPZ) ;

    //--------------------------------------------------------------------------
    // Declare vars
    //--------------------------------------------------------------------------

    int64_t p, j, i ;
    SPEX_matrix *b2 = NULL;   // b2 stores the solution of A*x
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);
    mpq_t scale; SPEX_MPQ_SET_NULL(scale);

    SPEX_CHECK (SPEX_mpq_init(temp));
    SPEX_CHECK(SPEX_mpq_init(scale));
    SPEX_CHECK (SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPQ, b->m, b->n,
        b->nzmax, false, true, option));


    //--------------------------------------------------------------------------
    // perform SPEX_mpq_addmul in loops to compute b2 = A'*x, where A' is the
    // scaled matrix with all entries in integer
    //--------------------------------------------------------------------------

    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            for (p = A->p[i]; p < A->p[i + 1]; p++)
            {
                // temp = A[p][i] (note this must be done seperately since A is
                // mpz and temp is mpq)
                SPEX_CHECK(SPEX_mpq_set_z(temp, A->x.mpz[p]));

                // temp = temp*x[i]
                SPEX_CHECK(SPEX_mpq_mul(temp, temp,
                                        SPEX_2D(x, i, j, mpq)));

                // b2[p] = b2[p]-temp
                SPEX_CHECK(SPEX_mpq_add(SPEX_2D(b2, A->i[p], j, mpq),
                                        SPEX_2D(b2, A->i[p], j, mpq),temp));
            }
        }
    }

    //--------------------------------------------------------------------------
    // Apply scales of A and b to b2 before comparing the b2 with scaled b'
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_mpq_div(scale, b->scale, A->scale));

    // Apply scaling factor, but ONLY if it is not 1
    int r;
    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, scale, 1, 1));
    if (r != 0)
    {
        for (i = 0; i < b2->m*b2->n; i++)
        {
            SPEX_CHECK(SPEX_mpq_mul(b2->x.mpq[i], b2->x.mpq[i], scale));
        }
    }

    //--------------------------------------------------------------------------
    // check if b==b2
    //--------------------------------------------------------------------------
    *Is_correct = true;
    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            // temp = b[i] (correct b)
            SPEX_CHECK(SPEX_mpq_set_z(temp, SPEX_2D(b, i, j, mpz)));

            // set check false if b!=b2
            SPEX_CHECK(SPEX_mpq_equal(&r, temp, SPEX_2D(b2, i, j, mpq)));
            if (r == 0)
            {
                *Is_correct = false;
                j = b->n;
                break;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Print info
    //--------------------------------------------------------------------------

    int pr = SPEX_OPTION_PRINT_LEVEL (option) ;
    if (*Is_correct)
    {
        SPEX_PR1 ("Solution is verified to be exact.\n") ;
    }
    else
    {
        // This can never happen.
        SPEX_PR1 ("ERROR! Solution is wrong. This is a bug; please "
                  "contact the authors of SPEX.\n") ;
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SPEX_FREE_ALL;
    return SPEX_OK;
}

