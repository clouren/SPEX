//------------------------------------------------------------------------------
// SPEX_Utilities/spex_gmp.h: definitions for SPEX_gmp.c
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// These macros are used by SPEX_gmp.c to create wrapper functions around all
// GMP functions used by SPEX, to safely handle out-of-memory conditions.
// They are placed in this separate #include file so that a future developer
// can use them to construct their own wrappers around GMP functions.  See
// SPEX_gmp.c for more details.

#ifndef SPEX_GMP_H
#define SPEX_GMP_H

#include "SPEX.h"

//------------------------------------------------------------------------------
//-------------------------functions for GMP wrapper----------------------------
//------------------------------------------------------------------------------

// uncomment this to print memory debugging info
// #define SPEX_GMP_MEMORY_DEBUG

#ifdef SPEX_GMP_MEMORY_DEBUG
void spex_gmp_dump ( void ) ;
#endif

#ifndef SPEX_GMP_LIST_INIT
// A size of 32 ensures that the list never needs to be increased in size.
// The test coverage suite in SPEX/Tcov reduces this initial size to
// exercise the code, in SPEX/Tcov/Makefile.
#define SPEX_GMP_LIST_INIT 32
#endif

void *spex_gmp_allocate (size_t size) ;

void spex_gmp_free (void *p, size_t size) ;

void *spex_gmp_reallocate (void *p_old, size_t old_size, size_t new_size );

void spex_gmp_failure (int status) ;

//------------------------------------------------------------------------------
// Field access macros for MPZ/MPQ/MPFR struct
//------------------------------------------------------------------------------
// FUTURE: make these accessible to the end user?

// (similar definition in gmp-impl.h and mpfr-impl.h)

#define SPEX_MPZ_SIZ(x)   ((x)->_mp_size)
#define SPEX_MPZ_PTR(x)   ((x)->_mp_d)
#define SPEX_MPZ_ALLOC(x) ((x)->_mp_alloc)
#define SPEX_MPQ_NUM(x)   mpq_numref(x)
#define SPEX_MPQ_DEN(x)   mpq_denref(x)
#define SPEX_MPFR_MANT(x) ((x)->_mpfr_d)
#define SPEX_MPFR_EXP(x)  ((x)->_mpfr_exp)
#define SPEX_MPFR_PREC(x) ((x)->_mpfr_prec)
#define SPEX_MPFR_SIGN(x) ((x)->_mpfr_sign)

/*re-define but same result: */
#define SPEX_MPFR_REAL_PTR(x) (&((x)->_mpfr_d[-1]))

/* Invalid exponent value (to track bugs...) */
#define SPEX_MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Macros to set the pointer in mpz_t/mpq_t/mpfr_t variable to NULL. It is best
 * practice to call these macros immediately after mpz_t/mpq_t/mpfr_t variable
 * is declared, and before the mp*_init function is called. It would help to
 * prevent error when SPEX_MP*_CLEAR is called before the variable is
 * successfully initialized.
 */

#define SPEX_MPZ_SET_NULL(x)                     \
{                                                \
    SPEX_MPZ_PTR(x) = NULL;                      \
    SPEX_MPZ_SIZ(x) = 0;                         \
    SPEX_MPZ_ALLOC(x) = 0;                       \
}

#define SPEX_MPQ_SET_NULL(x)                     \
{                                                \
    SPEX_MPZ_PTR(SPEX_MPQ_NUM(x)) = NULL;        \
    SPEX_MPZ_SIZ(SPEX_MPQ_NUM(x)) = 0;           \
    SPEX_MPZ_ALLOC(SPEX_MPQ_NUM(x)) = 0;         \
    SPEX_MPZ_PTR(SPEX_MPQ_DEN(x)) = NULL;        \
    SPEX_MPZ_SIZ(SPEX_MPQ_DEN(x)) = 0;           \
    SPEX_MPZ_ALLOC(SPEX_MPQ_DEN(x)) = 0;         \
}

#define SPEX_MPFR_SET_NULL(x)                    \
{                                                \
    SPEX_MPFR_MANT(x) = NULL;                    \
    SPEX_MPFR_PREC(x) = 0;                       \
    SPEX_MPFR_SIGN(x) = 1;                       \
    SPEX_MPFR_EXP(x) = SPEX_MPFR_EXP_INVALID;    \
}

/* GMP does not give a mechanism to tell a user when an mpz, mpq, or mpfr
 * item has been cleared; thus, if mp*_clear is called on an object that
 * has already been cleared, gmp will crash. It is also not possible to
 * set a mp*_t = NULL. Thus, this mechanism modifies the internal GMP
 * size of entries to avoid crashing in the case that a mp*_t is cleared
 * multiple times.
 */

#define SPEX_MPZ_CLEAR(x)                        \
{                                                \
    if ((x) != NULL && SPEX_MPZ_PTR(x) != NULL)  \
    {                                            \
        mpz_clear(x);                            \
        SPEX_MPZ_SET_NULL(x);                    \
    }                                            \
}

#define SPEX_MPQ_CLEAR(x)                        \
{                                                \
    SPEX_MPZ_CLEAR(SPEX_MPQ_NUM(x));             \
    SPEX_MPZ_CLEAR(SPEX_MPQ_DEN(x));             \
}

#define SPEX_MPFR_CLEAR(x)                       \
{                                                \
    if ((x) != NULL && SPEX_MPFR_MANT(x) != NULL)\
    {                                            \
        mpfr_clear(x);                           \
        SPEX_MPFR_SET_NULL(x);                   \
    }                                            \
}

//------------------------------------------------------------------------------
// GMP/MPFR wrapper macros
//------------------------------------------------------------------------------

#define SPEX_GMP_WRAPPER_START                                          \
{                                                                       \
    spex_gmp_nmalloc = 0 ;                                              \
    /* setjmp returns 0 if called from here, or > 0 if from longjmp */  \
    int spex_gmp_status = setjmp (spex_gmp_environment) ;               \
    if (spex_gmp_status != 0)                                           \
    {                                                                   \
        /* failure from longjmp */                                      \
        spex_gmp_failure (spex_gmp_status) ;                            \
        return (SPEX_OUT_OF_MEMORY) ;                                   \
    }                                                                   \
}

#define SPEX_GMPZ_WRAPPER_START(x)                                      \
{                                                                       \
    spex_gmpz_archive = (mpz_t *) x;                                    \
    spex_gmpz_archive2 = NULL;                                          \
    spex_gmpq_archive = NULL;                                           \
    spex_gmpfr_archive = NULL;                                          \
    SPEX_GMP_WRAPPER_START;                                             \
}

#define SPEX_GMPZ_WRAPPER_START2(x,y)                                   \
{                                                                       \
    spex_gmpz_archive  = (mpz_t *) x;                                   \
    spex_gmpz_archive2 = (mpz_t *) y;                                   \
    spex_gmpq_archive = NULL;                                           \
    spex_gmpfr_archive = NULL;                                          \
    SPEX_GMP_WRAPPER_START;                                             \
}

#define SPEX_GMPQ_WRAPPER_START(x)                                      \
{                                                                       \
    spex_gmpz_archive = NULL;                                           \
    spex_gmpz_archive2 = NULL;                                          \
    spex_gmpq_archive =(mpq_t *) x;                                     \
    spex_gmpfr_archive = NULL;                                          \
    SPEX_GMP_WRAPPER_START;                                             \
}

#define SPEX_GMPFR_WRAPPER_START(x)                                     \
{                                                                       \
    spex_gmpz_archive = NULL;                                           \
    spex_gmpz_archive2 = NULL;                                          \
    spex_gmpq_archive = NULL;                                           \
    spex_gmpfr_archive = (mpfr_t *) x;                                  \
    SPEX_GMP_WRAPPER_START;                                             \
}

#define SPEX_GMP_WRAPPER_FINISH                                         \
{                                                                       \
    /* clear (but do not free) the list.  The caller must ensure */     \
    /* the result is eventually freed. */                               \
    spex_gmpz_archive = NULL ;                                          \
    spex_gmpz_archive2 = NULL;                                          \
    spex_gmpq_archive = NULL ;                                          \
    spex_gmpfr_archive = NULL ;                                         \
    spex_gmp_nmalloc = 0 ;                                              \
}

// free a block of memory, and also remove it from the archive if it's there
#define SPEX_GMP_SAFE_FREE(p)                                           \
{                                                                       \
    if (spex_gmpz_archive != NULL)                                      \
    {                                                                   \
        if (p == SPEX_MPZ_PTR(*spex_gmpz_archive))                      \
        {                                                               \
            SPEX_MPZ_PTR(*spex_gmpz_archive) = NULL ;                   \
        }                                                               \
    }                                                                   \
    else if (spex_gmpz_archive2 != NULL)                                \
    {                                                                   \
        if (p == SPEX_MPZ_PTR(*spex_gmpz_archive2))                     \
        {                                                               \
            SPEX_MPZ_PTR(*spex_gmpz_archive2) = NULL ;                  \
        }                                                               \
    }                                                                   \
    else if (spex_gmpq_archive != NULL)                                 \
    {                                                                   \
        if (p == SPEX_MPZ_PTR(SPEX_MPQ_NUM(*spex_gmpq_archive)))        \
        {                                                               \
            SPEX_MPZ_PTR(SPEX_MPQ_NUM(*spex_gmpq_archive)) = NULL ;     \
        }                                                               \
        if (p == SPEX_MPZ_PTR(SPEX_MPQ_DEN(*spex_gmpq_archive)))        \
        {                                                               \
            SPEX_MPZ_PTR(SPEX_MPQ_DEN(*spex_gmpq_archive)) = NULL ;     \
        }                                                               \
    }                                                                   \
    else if (spex_gmpfr_archive != NULL)                                \
    {                                                                   \
        if (p == SPEX_MPFR_REAL_PTR(*spex_gmpfr_archive))               \
        {                                                               \
            SPEX_MPFR_MANT(*spex_gmpfr_archive) = NULL ;                \
        }                                                               \
    }                                                                   \
    SPEX_FREE (p) ;                                                     \
}


#endif
