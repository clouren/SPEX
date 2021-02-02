//------------------------------------------------------------------------------
// SPEX_Left_LU/Source/spex_left_lu_internal: include file for internal use in
// SPEX_Left_LU
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_Left_LU.h instead.

#ifndef SPEX_LEFT_LU_INTERNAL_H
#define SPEX_LEFT_LU_INTERNAL_H

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-value"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------C Libraries------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Standard C libraries
#include <setjmp.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
// debugging
//------------------------------------------------------------------------------

#ifdef SPEX_DEBUG

#ifdef MATLAB_MEX_FILE

#define ASSERT(x)                                                             \
{                                                                             \
    if (!(x))                                                                 \
    {                                                                         \
        mexErrMsgTxt ("assertion failed: %s line %d\n", __FILE__, __LINE__) ; \
    }                                                                         \
}

#else

#include <assert.h>
#define ASSERT(x) assert (x)

#endif

#else

// debugging disabled
#define ASSERT(x)

#endif

//------------------------------------------------------------------------------
//-------------------------Common Macros----------------------------------------
//------------------------------------------------------------------------------

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

#define SPEX_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SPEX_FLIP(i) (-(i)-2)
#define SPEX_UNFLIP(i) (((i) < 0) ? SPEX_FLIP(i) : (i))
#define SPEX_MARKED(Ap,j) (Ap [j] < 0)
#define SPEX_MARK(Ap,j) { Ap [j] = SPEX_FLIP (Ap [j]) ; }

// SPEX_CHECK(method) is a macro that calls a SPEX LU method and checks the
// status; if a failure occurs, it frees all allocated workspace and returns
// the error status to the caller.  To use SPEX_CHECK, the #include'ing file
// must declare a SPEX_info info, and must define SPEX_FREE_ALL as a macro that
// frees all workspace if an error occurs. The method can be a scalar as well,
// so that SPEX_CHECK(info) works.

// the default is to free nothing
#ifndef SPEX_FREE_ALL
#define SPEX_FREE_ALL
#endif

#define SPEX_CHECK(method)      \
{                               \
    info = (method) ;           \
    if (info != SPEX_OK)        \
    {                           \
        SPEX_FREE_ALL ;         \
        return (info) ;         \
    }                           \
}

#include "SPEX_Left_LU.h"

//------------------------------------------------------------------------------
// printing control
//------------------------------------------------------------------------------

// SPEX_LU uses SuiteSparse_config.printf_func instead of a mere call to printf
// (the default function is printf, or mexPrintf when in MATLAB).  If this
// function pointer is NULL, no printing is done.

#define SPEX_PRINTF(...)                                    \
{                                                           \
    if (SuiteSparse_config.printf_func != NULL)             \
    {                                                       \
        SuiteSparse_config.printf_func (__VA_ARGS__) ;      \
    }                                                       \
}

#define SPEX_PR1(...) { if (pr >= 1) SPEX_PRINTF (__VA_ARGS__) }
#define SPEX_PR2(...) { if (pr >= 2) SPEX_PRINTF (__VA_ARGS__) }
#define SPEX_PR3(...) { if (pr >= 3) SPEX_PRINTF (__VA_ARGS__) }

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------functions for GMP wrapper----------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// uncomment this to print memory debugging info
// #define SPEX_GMP_MEMORY_DEBUG


#define SPEX_DEFAULT_TOL 1

// Check parameter. If this = 1 then the solution to the system is checked
// for accuracy
#define SPEX_DEFAULT_CHECK false

// Pivoting scheme used for SPEX LU.
//  SPEX_SMALLEST = 0,              Smallest pivot
//  SPEX_DIAGONAL = 1,              Diagonal pivoting
//  SPEX_FIRST_NONZERO = 2,         First nonzero per column chosen as pivot
//  SPEX_TOL_SMALLEST = 3,          Diagonal pivoting with tolerance for small
//  SPEX_TOL_LARGEST = 4,           Diagonal pivoting with tolerance for large
//  SPEX_LARGEST = 5                Largest pivot
#define SPEX_DEFAULT_PIVOT SPEX_TOL_SMALLEST

// Column ordering used.
//  SPEX_NO_ORDERING = 0,           None: Not recommended for sparse matrices
//  SPEX_COLAMD = 1,                COLAMD: Default
//  SPEX_AMD = 2                    AMD
#define SPEX_DEFAULT_ORDER SPEX_COLAMD

// Defines printing to be done
#define SPEX_DEFAULT_PRINT_LEVEL 0

// MPFR precision used (quad is default)
#define SPEX_DEFAULT_PRECISION 128

//------------------------------------------------------------------------------
// Type of MPFR rounding used.
//------------------------------------------------------------------------------

// The MPFR library utilizes an internal rounding scheme. The options are
//  MPFR_RNDN: round to nearest (roundTiesToEven in IEEE 754-2008),
//  MPFR_RNDZ: round toward zero (roundTowardZero in IEEE 754-2008),
//  MPFR_RNDU: round toward plus infinity (roundTowardPositive in
//             IEEE 754-2008),
//  MPFR_RNDD: round toward minus infinity (roundTowardNegative in
//             IEEE 754-2008),
//  MPFR_RNDA: round away from zero.
//  MPFR_RNDF: faithful rounding. This is not stable
//
// SPEX LU utilizes MPFR_RNDN by default.

#define SPEX_DEFAULT_MPFR_ROUND MPFR_RNDN

//------------------------------------------------------------------------------
// Macros to utilize the default if option is NULL
//------------------------------------------------------------------------------

#define SPEX_OPTION(option,parameter,default_value) \
    ((option == NULL) ? (default_value) : (option->parameter))

#define SPEX_OPTION_TOL(option) \
    SPEX_OPTION (option, tol, SPEX_DEFAULT_TOL)

#define SPEX_OPTION_CHECK(option) \
    SPEX_OPTION (option, check, false)

#define SPEX_OPTION_PIVOT(option) \
    SPEX_OPTION (option, pivot, SPEX_DEFAULT_PIVOT)

#define SPEX_OPTION_ORDER(option) \
    SPEX_OPTION (option, order, SPEX_DEFAULT_ORDER)

#define SPEX_OPTION_PREC(option) \
    SPEX_OPTION (option, prec, SPEX_DEFAULT_PRECISION)

#define SPEX_OPTION_PRINT_LEVEL(option) \
    SPEX_OPTION (option, print_level, SPEX_DEFAULT_PRINT_LEVEL)

#define SPEX_OPTION_ROUND(option) \
    SPEX_OPTION (option, round, SPEX_DEFAULT_MPFR_ROUND)

//------------------------------------------------------------------------------
// Field access macros for MPZ/MPQ/MPFR struct
//------------------------------------------------------------------------------

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

#define SPEX_MPZ_SET_NULL(x)                \
    SPEX_MPZ_PTR(x) = NULL;                 \
    SPEX_MPZ_SIZ(x) = 0;                    \
    SPEX_MPZ_ALLOC(x) = 0;

#define SPEX_MPQ_SET_NULL(x)                     \
    SPEX_MPZ_PTR(SPEX_MPQ_NUM(x)) = NULL;        \
    SPEX_MPZ_SIZ(SPEX_MPQ_NUM(x)) = 0;           \
    SPEX_MPZ_ALLOC(SPEX_MPQ_NUM(x)) = 0;         \
    SPEX_MPZ_PTR(SPEX_MPQ_DEN(x)) = NULL;        \
    SPEX_MPZ_SIZ(SPEX_MPQ_DEN(x)) = 0;           \
    SPEX_MPZ_ALLOC(SPEX_MPQ_DEN(x)) = 0;

#define SPEX_MPFR_SET_NULL(x)               \
    SPEX_MPFR_MANT(x) = NULL;               \
    SPEX_MPFR_PREC(x) = 0;                  \
    SPEX_MPFR_SIGN(x) = 1;                  \
    SPEX_MPFR_EXP(x) = SPEX_MPFR_EXP_INVALID;

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

#define SPEX_MPQ_CLEAR(x)                   \
{                                           \
    SPEX_MPZ_CLEAR(SPEX_MPQ_NUM(x));        \
    SPEX_MPZ_CLEAR(SPEX_MPQ_DEN(x));        \
}

#define SPEX_MPFR_CLEAR(x)                        \
{                                                 \
    if ((x) != NULL && SPEX_MPFR_MANT(x) != NULL) \
    {                                             \
        mpfr_clear(x);                            \
        SPEX_MPFR_SET_NULL(x);                    \
    }                                             \
}


// ============================================================================
//                           Internal Functions
// ============================================================================


/* Purpose: This function performs sparse REF forward substitution. This is
 * essentially the same as the sparse REF triangular solve applied to each
 * column of the right hand side vectors. Like the normal one, this function
 * expects that the matrix x is dense. As a result,the nonzero pattern is not
 * computed and each nonzero in x is iterated across.  The system to solve is
 * LDx = x.  On output, the mpz_t** x structure is modified.
 */
SPEX_info spex_left_lu_forward_sub
(
    const SPEX_matrix *L,   // lower triangular matrix
    SPEX_matrix *x,         // right hand side matrix
    const SPEX_matrix *rhos // sequence of pivots used in factorization
);

/* Purpose: This function performs sparse REF backward substitution. In essense
 * it solves the sysem Ux = x. Note that prior to this, we expect x to be
 * multiplied by the determinant of A.  The input argument bx is modified on
 * output.
 */
SPEX_info spex_left_lu_back_sub  // performs sparse REF backward substitution
(
    const SPEX_matrix *U,   // input upper triangular matrix
    SPEX_matrix *bx        // right hand side matrix of size n*numRHS
)  ;


/* Purpose: This function performs a depth first search of the graph of the
 * matrix starting at node j. The output of this function is the set of nonzero
 * indices in the xi vector.
 */
void spex_left_lu_dfs // performs a dfs of the graph of the matrix starting at node j
(
    int64_t *top,          // beginning of stack
    int64_t j,             // What node to start DFS at
    SPEX_matrix* L,        // matrix which represents the Graph of L
    int64_t* xi,           // the nonzero pattern
    int64_t* pstack,       // workspace vector
    const int64_t* pinv   // row permutation
);


/* This function performs the pivoting for the SPEX LU factorization.
 * The optional Order is:
 *     0: Smallest pivot
 *     1: Natural/Diagonal pivoting
 *     2: Choose first nonzero (not recommended, for comparison only)
 *     3: Diagonal with tolerance and smallest pivot (default)
 *     4: Diagonal with tolerance and largest pivoting
 *     5: Largest pivot (not recommended, for comparison only)
 *
 * On output, the pivs, pinv, and row_perm arrays and rhos matrix are all modified.
 */
SPEX_info spex_left_lu_get_pivot
(
    int64_t *pivot,      // found index of pivot entry
    SPEX_matrix* x,      // kth column of L and U
    int64_t* pivs,       // vector indicating which rows have been pivotal
    int64_t n,           // dimension of the problem
    int64_t top,         // nonzero pattern is located in xi[top..n-1]
    int64_t* xi,         // nonzero pattern of x
    int64_t col,         // current column of A (real kth column i.e., q[k])
    int64_t k,           // iteration of the algorithm
    SPEX_matrix* rhos,   // vector of pivots
    int64_t* pinv,       // row permutation
    int64_t* row_perm,   // opposite of pinv. if pinv[i] = j then row_perm[j] = i
    const SPEX_options *option // command option
);

/* Purpose: This function selects the pivot element as the largest in the
 * column. This is activated if the user sets option->pivot = SPEX_LARGEST.
 * NOTE: This pivoting scheme is NOT recommended for SPEX LU.  On output, the
 * index of the largest pivot is returned.
 */
SPEX_info spex_left_lu_get_largest_pivot
(
    int64_t *pivot,         // the index of largest pivot
    SPEX_matrix* x,         // kth column of L and U
    int64_t* pivs,          // vector which indicates whether each row has been pivotal
    int64_t n,              // dimension of problem
    int64_t top,            // nonzero pattern is located in xi[top..n-1]
    int64_t* xi             // nonzero pattern of x
);

/* This function obtains the first eligible nonzero pivot.  This is enabled if
 * the user sets option->pivot = SPEX_FIRST_NONZERO.  NOTE: This pivoting
 * scheme is not recommended.  On output, the kth pivot is returned.
 */
SPEX_info spex_left_lu_get_nonzero_pivot // find the first eligible nonzero pivot
(
    int64_t *pivot,      // the index of first eligible nonzero pivot
    SPEX_matrix* x,      // kth column of L and U
    int64_t* pivs,       // vector indicating which rows are pivotal
    int64_t n,           // size of x
    int64_t top,         // nonzero pattern is located in xi[top..n-1]
    int64_t* xi          // nonzero pattern of x
);

/* Purpose: This function selects the pivot element as the smallest in the
 * column. This is activated by default or if the user sets option->pivot =
 * SPEX_SMALLEST.  NOTE: This is the recommended pivoting scheme for SPEX LU.
 * On output, the index of kth pivot is returned.
 */
SPEX_info spex_left_lu_get_smallest_pivot
(
    int64_t *pivot,         // index of smallest pivot
    SPEX_matrix *x,         // kth column of L and U
    int64_t* pivs,          // vector indicating whether each row has been pivotal
    int64_t n,              // dimension of problem
    int64_t top,            // nonzeros are stored in xi[top..n-1]
    int64_t* xi             // nonzero pattern of x
);


/* Purpose: This function permutes b for forward substitution.
 * That is, b = P'*b.
 */
SPEX_info spex_left_lu_permute_b
(
    SPEX_matrix **b_handle,     // permuted RHS vector
    const SPEX_matrix *b2,      // unpermuted RHS vector (not modified)
    const int64_t *pinv,        // inverse row permutation
    const SPEX_options* option
);

/* Purpose: SPEX_permute_x permutes x to get it back in its original form.
 * That is x = Q*x.
 */
SPEX_info spex_left_lu_permute_x
(
    SPEX_matrix **x_handle,    // permuted Solution vector
    SPEX_matrix *x2,           // unpermuted Solution vector (not modified)
    SPEX_LU_analysis *S,  // symbolic analysis with the column ordering Q
    const SPEX_options* option  // Command options
                          // has been checked in the only caller SPEX_LU_solve
) ;

/* Purpose: This function computes the reach of column k of A on the graph of L
 * mathematically that is: xi = Reach(A(:,k))_G_L.
 */
void spex_left_lu_reach    // compute the reach of column k of A on the graph of L
(
    int64_t *top,
    SPEX_matrix* L,         // matrix representing graph of L
    const SPEX_matrix* A,   // input matrix
    int64_t k,              // column of A of interest
    int64_t* xi,            // nonzero pattern
    const int64_t* pinv     // row permutation
)  ;

/* Purpose: This function performs the sparse REF triangular solve; that is,
 * (LD) x = A(:,k). The algorithm is described in the paper; however in essence
 * it computes the nonzero pattern xi, then performs a sequence of IPGE
 * operations on the nonzeros to obtain their final value. All operations are
 * gauranteed to be integral. There are various enhancements in this code used
 * to reduce the overall cost of the operations and minimize operations as much
 * as possible.
 */
SPEX_info spex_left_lu_ref_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,      // Output the beginning of nonzero pattern
    SPEX_matrix* L,           // partial L matrix
    const SPEX_matrix* A,     // input matrix
    int64_t k,                // iteration of algorithm
    int64_t* xi,              // nonzero pattern vector
    const int64_t* q,         // column permutation
    SPEX_matrix* rhos,        // sequence of pivots
    const int64_t* pinv,      // inverse row permutation
    const int64_t* row_perm,  // row permutation
    int64_t* h,               // history vector
    SPEX_matrix* x            // solution of system ==> kth column of L and U
);


// (void *) pointer to the values of A.  A must be non-NULL with a valid type
#define SPEX_X(A)                                                           \
    ((A->type == SPEX_MPZ  ) ? (void *) A->x.mpz   :                        \
    ((A->type == SPEX_MPQ  ) ? (void *) A->x.mpq   :                        \
    ((A->type == SPEX_MPFR ) ? (void *) A->x.mpfr  :                        \
    ((A->type == SPEX_INT64) ? (void *) A->x.int64 : (void *) A->x.fp64))))


// return an error if A->kind (csc, triplet, dense) is wrong
#define SPEX_REQUIRE_KIND(A,required_kind) \
    if (A == NULL || A->kind != required_kind) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_KIND(A,required_kind) \
    ASSERT (A != NULL && A->kind == required_kind)

// return an error if A->type (mpz, mpq, mpfr, int64, or double) is wrong
#define SPEX_REQUIRE_TYPE(A,required_type) \
    if (A == NULL || A->type != required_type) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_TYPE(A,required_type) \
    ASSERT (A != NULL && A->type == required_type)

// return an error if A->kind or A->type is wrong
#define SPEX_REQUIRE(A,required_kind,required_type)     \
    SPEX_REQUIRE_KIND (A,required_kind) ;               \
    SPEX_REQUIRE_TYPE (A,required_type) ;

#define ASSERT_MATRIX(A,required_kind,required_type)    \
    ASSERT_KIND (A,required_kind) ;                     \
    ASSERT_TYPE (A,required_type) ;

#endif

