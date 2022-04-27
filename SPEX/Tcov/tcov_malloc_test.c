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
    const SPEX_matrix *A,         // Input matrix
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
                info = SPEX_INCORRECT;
                j = b->n;
                break;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Print info
    //--------------------------------------------------------------------------

    int pr = SPEX_OPTION_PRINT_LEVEL (option) ;
    if (info == SPEX_OK)
    {
        SPEX_PR1 ("Solution is verified to be exact.\n") ;
    }
    else if (info == SPEX_INCORRECT)
    {
        // This can never happen.
        SPEX_PR1 ("ERROR! Solution is wrong. This is a bug; please "
                  "contact the authors of SPEX.\n") ;
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SPEX_FREE_ALL;
    return info;
}

//------------------------------------------------------------------------------
// spex_update_verify.c: verify if A=LD^(-1)U after factorization update
//------------------------------------------------------------------------------

/* Purpose: This function is to verify if A=L(P,:)D^(-1)U(:,Q) after
 * factorization update. This is done by solving LD^(-1)U*x=b via the updated
 * factorization, and check if A*x=b holds rational-arthmetically. This
 * function is provided here only used for debugging purposes, as the routines
 * within SPEX are gauranteed to be exact.
 */


#undef SPEX_FREE_ALL
#define SPEX_FREE_ALL                 \
    SPEX_matrix_free(&b, option);     \
    SPEX_matrix_free(&x, option);     \
    SPEX_matrix_free(&b2, option);    \
    SPEX_MPQ_CLEAR(temp);


SPEX_info spex_update_verify
(
    SPEX_factorization *F,// LU factorization of A
    const SPEX_matrix *A,     // Input matrix
    const SPEX_options *option// command options
)
{
    SPEX_info info;
    int64_t tmp, i, n = F->L->n;
    int r;
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);
    SPEX_matrix *b = NULL; // the dense right-hand-side matrix to be generated
    SPEX_matrix *x = NULL; // the dense solution matrix to be generated
    SPEX_matrix *b2 = NULL; // the dense matrix to store the result of A*x

    SPEX_CHECK(SPEX_mpq_init(temp));
    SPEX_CHECK(SPEX_matrix_allocate(&b , SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option));
    SPEX_CHECK(SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPQ, n, 1, n, false,
        true, option));

    // -------------------------------------------------------------------------
    // generate random right-hand-size vector
    // -------------------------------------------------------------------------
    // initialize random number generator
    int seed = 10;
    srand(seed);
    for (i = 0; i < n; i++)
    {
        tmp = i+1;//rand(); //TODO
        SPEX_CHECK(SPEX_mpz_set_si(b->x.mpz[i], tmp));
    }

    // -------------------------------------------------------------------------
    // solve LD^(-1)Ux = b for x
    // -------------------------------------------------------------------------
    SPEX_CHECK(SPEX_Update_solve(&x, F, b, option));

    // -------------------------------------------------------------------------
    // compute b2 = A*x
    // -------------------------------------------------------------------------
    for (i = 0; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpz_sgn(&r, x->x.mpz[i]));
        if (r == 0) { continue;}

        for (int64_t p = 0; p < A->v[i]->nz; p++)
        {
            int64_t j = A->v[i]->i[p];
            SPEX_CHECK(SPEX_mpq_set_z(temp, A->v[i]->x[p]));
            // b2[j] += x[i]*A(j,i)
            SPEX_CHECK(SPEX_mpq_mul(temp, temp, x->x.mpq[i]));
            SPEX_CHECK(SPEX_mpq_add(b2->x.mpq[j], b2->x.mpq[j], temp));
        }
    }
    //--------------------------------------------------------------------------
    // Apply scales of A and b to b2 before comparing the b2 with scaled b'
    //--------------------------------------------------------------------------
    SPEX_CHECK(SPEX_mpq_div(temp, b->scale, A->scale));

    // Apply scaling factor, but ONLY if it is not 1
    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, temp, 1, 1));
    if (r != 0)
    {
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpq_mul(b2->x.mpq[i], b2->x.mpq[i], temp));
        }
    }

    // -------------------------------------------------------------------------
    // check if b2 == b
    // -------------------------------------------------------------------------
    for (i = 0; i < n; i++)
    {
        // temp = b[i] (correct b)
        SPEX_CHECK(SPEX_mpq_set_z(temp, b->x.mpz[i]));

        // set check false if b!=b2
        SPEX_CHECK(SPEX_mpq_equal(&r, temp, b2->x.mpq[i]));
        if (r == 0)
        {
            info = SPEX_INCORRECT;
            break;
        }
    }

    //--------------------------------------------------------------------------
    // Print info
    //--------------------------------------------------------------------------

    int pr = SPEX_OPTION_PRINT_LEVEL (option) ;
    if (info == SPEX_OK)
    {
        SPEX_PR1 ("Factorization is verified to be correct and exact.\n") ;
    }
    else if (info == SPEX_INCORRECT)
    {
        // This can never happen.
        SPEX_PR1 ("ERROR! Factorization is wrong. This is a bug; please "
                  "contact the authors of SPEX.\n") ;
    }

    SPEX_FREE_ALL;
    return info;
}
