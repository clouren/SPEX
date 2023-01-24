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
// spex_gmp environment
//------------------------------------------------------------------------------

#include <setjmp.h>

typedef struct
{
    jmp_buf environment ;   // for setjmp and longjmp
    int64_t nmalloc ;       // # of malloc'd objects in spex_gmp->list
    int64_t nlist ;         // size of the spex_gmp->list
    void **list ;           // list of malloc'd objects
    mpz_t  *mpz_archive  ;  // current mpz object
    mpz_t  *mpz_archive2 ;  // current second mpz object
    mpq_t  *mpq_archive  ;  // current mpq object
    mpfr_t *mpfr_archive ;  // current mpfr object
}
spex_gmp_t ;

#ifndef SPEX_GMP_LIST_INIT
// Initial size of the spex_gmp->list.  A size of 32 ensures that the list
// never needs to be increased in size (at least in practice; it is possible
// that GMP or MPFR could exceed this size).  The test coverage suite in
// SPEX/Tcov reduces this initial size to exercise the code, in
// SPEX/Tcov/Makefile.
#define SPEX_GMP_LIST_INIT 32
#endif

// for debugging only:
SUITESPARSE_PUBLIC int64_t spex_gmp_ntrials ;

//------------------------------------------------------------------------------
// SPEX GMP functions
//------------------------------------------------------------------------------

// uncomment this to print memory debugging info
// #define SPEX_GMP_MEMORY_DEBUG

int spex_gmp_initialize (void) ;

void spex_gmp_finalize (void) ;

spex_gmp_t *spex_gmp_get (void) ;

void *spex_gmp_allocate (size_t size) ;

void spex_gmp_free (void *p, size_t size) ;

void *spex_gmp_reallocate (void *p_old, size_t old_size, size_t new_size );

#ifdef SPEX_GMP_MEMORY_DEBUG
void spex_gmp_dump ( void ) ;
#endif

SPEX_info spex_gmp_failure (int status) ;

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

#define SPEX_GMP_WRAPPER_START_HELPER(z1,z2,q,fr)                       \
    spex_gmp_t *spex_gmp = spex_gmp_get ( ) ;                           \
    spex_gmp->mpz_archive  = (mpz_t  *) z1 ;                            \
    spex_gmp->mpz_archive2 = (mpz_t  *) z2 ;                            \
    spex_gmp->mpq_archive  = (mpq_t  *) q  ;                            \
    spex_gmp->mpfr_archive = (mpfr_t *) fr ;                            \
    /* setjmp returns 0 if called from here, or > 0 if from longjmp */  \
    int status = setjmp (spex_gmp->environment) ;                       \
    if (status != 0)                                                    \
    {                                                                   \
        /* failure from longjmp */                                      \
        return (spex_gmp_failure (status)) ;                            \
    }

#define SPEX_GMP_WRAPPER_START                                          \
    SPEX_GMP_WRAPPER_START_HELPER (NULL, NULL, NULL, NULL) ;

#define SPEX_GMPZ_WRAPPER_START(z1)                                     \
    SPEX_GMP_WRAPPER_START_HELPER (z1, NULL, NULL, NULL) ;

#define SPEX_GMPZ_WRAPPER_START2(z1,z2)                                 \
    SPEX_GMP_WRAPPER_START_HELPER (z1, z2, NULL, NULL) ;

#define SPEX_GMPQ_WRAPPER_START(q)                                      \
    SPEX_GMP_WRAPPER_START_HELPER (NULL, NULL, q, NULL) ;

#define SPEX_GMPFR_WRAPPER_START(fr)                                    \
    SPEX_GMP_WRAPPER_START_HELPER (NULL, NULL, NULL, fr) ;

#define SPEX_GMP_WRAPPER_FINISH                                         \
    spex_gmp->nmalloc = 0 ;                                             \
    spex_gmp->mpz_archive  = NULL ;                                     \
    spex_gmp->mpz_archive2 = NULL ;                                     \
    spex_gmp->mpq_archive  = NULL ;                                     \
    spex_gmp->mpfr_archive = NULL ;

// free a block of memory, and also remove it from the archive if it's there
#define SPEX_GMP_SAFE_FREE(p)                                           \
{                                                                       \
    if (spex_gmp->mpz_archive != NULL)                                  \
    {                                                                   \
        if (p == SPEX_MPZ_PTR(*(spex_gmp->mpz_archive)))                \
        {                                                               \
            SPEX_MPZ_PTR(*(spex_gmp->mpz_archive)) = NULL ;             \
        }                                                               \
    }                                                                   \
    else if (spex_gmp->mpz_archive2 != NULL)                            \
    {                                                                   \
        if (p == SPEX_MPZ_PTR(*spex_gmp->mpz_archive2))                 \
        {                                                               \
            SPEX_MPZ_PTR(*spex_gmp->mpz_archive2) = NULL ;              \
        }                                                               \
    }                                                                   \
    else if (spex_gmp->mpq_archive != NULL)                             \
    {                                                                   \
        if (p == SPEX_MPZ_PTR(SPEX_MPQ_NUM(*spex_gmp->mpq_archive)))    \
        {                                                               \
            SPEX_MPZ_PTR(SPEX_MPQ_NUM(*spex_gmp->mpq_archive)) = NULL ; \
        }                                                               \
        if (p == SPEX_MPZ_PTR(SPEX_MPQ_DEN(*spex_gmp->mpq_archive)))    \
        {                                                               \
            SPEX_MPZ_PTR(SPEX_MPQ_DEN(*spex_gmp->mpq_archive)) = NULL ; \
        }                                                               \
    }                                                                   \
    else if (spex_gmp->mpfr_archive != NULL)                            \
    {                                                                   \
        if (p == SPEX_MPFR_REAL_PTR(*spex_gmp->mpfr_archive))           \
        {                                                               \
            SPEX_MPFR_MANT(*spex_gmp->mpfr_archive) = NULL ;            \
        }                                                               \
    }                                                                   \
    SPEX_FREE (p) ;                                                     \
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------------GMP/MPFR wrapper macros------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// GMP/MPFR wrapper macros that already incorporate the SPEX_CHECK macro, which   
// is use to prevent segmentation faults in case gmp runs out of memory inside 
// a gmp call.
//------------------------------------------------------------------------------


#define SPEX_MPZ_MUL(a,b,c) SPEX_CHECK( SPEX_mpz_mul(a, b, c) )
#define SPEX_MPZ_MUL_SI(a,b,c) SPEX_CHECK( SPEX_mpz_mul_si(a, b, c) )


#endif
