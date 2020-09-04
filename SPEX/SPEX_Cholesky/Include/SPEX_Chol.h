//------------------------------------------------------------------------------
// SPEX_Cholesky/Include/SPEX_Cholesky.h: user #include file for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_CHOLESKY_H
#define SPEX_CHOLESKY_H


// This software package exactly solves a sparse symmetric positive definite (SPD)
// system of linear equations using one of two Integer-Preserving Cholesky factorizations. 
// This code accompanies the paper (to be submitted to ACM TOMs)

//    "Algorithm 1XXX: Exactly Solving Sparse Symmetric Positive Definite Linear Systems
//     via SPEX Cholesky factorization," C. Lourenco, E. Moreno-Centeno, T. Davis,
//     to be submitted ACM TOMS.

//   The theory associated with this paper is found at:

//    "A Sparse Integer-Preserving Cholesky Factorization Computed in Time Proportional
//     to Arithmetic Work", C. Lourenco, E. Moreno-Centeno, under submission, SIMAX.

//    If you use this code, you must first download and install the GMP, 
//    MPFR, SPEX_Left_LU, AMD, and COLAMD libraries. 
//   
//   GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/
//
//   SPEX_Left_LU, AMD, and COLAMD are distributed along with SPEX_Cholesky. The easiest
//   way ensure these dependencies are met is to only access this package through 
//   the SPEX repository.
//
//   All of these codes are components of the SPEX software library. This code may
//   be found at:
//              https://github.com/clouren/spex
//              www.suitesparse.com
//  
//


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Erick Moreno-Centeno, Timothy Davis, Jinhao Chen

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Chris Lourenco (chrisjlourenco@gmail.com, lourenco@usna.edu)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX_Cholesky is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Christopher Lourenco, Erick
// Moreno-Centeno, Timothy A. Davis and Jinhao Chen. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX Cholesky is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package solves the SPD linear system Ax = b exactly. The key
//    property of this package is that it can exactly solve ALL SPD input systems. 
//    The input matrix and right hand side vectors are stored as either integers, double
//    precision numbers, multiple precision floating points (through the mpfr
//    library) or as rational numbers (as a collection of numerators and
//    denominators using the GMP mpq_t data structure). Appropriate routines
//    within the code transform the input into an integral matrix in compressed
//    column form.

//    This package computes the factorization PAP' = LDL'. Note that we store the
//    "functional" form of the factorization by only storing the matrix L. The user
//    is given some freedom to select the permutation matrix P. The
//    recommended default settings select P using the AMD ordering.
//    Alternative strategies allowed to select P include the COLAMD 
//    ordering or no column permutation (P=I).  

//    The factor L is computed via integer preserving operations via
//    integer-preserving Gaussian elimination. The key part of this algorithm
//    is a REF Sparse triangular solve function which exploits sparsity and symmetry to
//    reduce the number of operations that must be performed.

//    Once L is computed, a simplified version of the triangular solve
//    is performed which assumes the vector b is dense. The final solution
//    vector x is gauranteed to be exact. This vector can be output in one of
//    three ways: 1) full precision rational arithmetic (as a sequence of
//    numerators and denominators) using the GMP mpq_t data type, 2) double
//    precision while not exact will produce a solution accurate to machine
//    roundoff unless the size of the associated solution exceeds double
//    precision (i.e., the solution is 10^500 or something), 3) variable
//    precision floating point using the GMP mpfr_t data type. The associated
//    precision is user defined.


// Here are the things to do
// TODO: Replace all prints with SPEX_print
// TODO: Documentation
// TODO: Tcov
// TODO: Change outputs to SPEX_info
// TODO: Make internal functions and SPEX_Chol_internal
// TODO: Put const where appropriate

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX Cholesly------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SPEX_Util.h"

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_CHOL_VERSION "1.0.0"
#define SPEX_CHOL_VERSION_MAJOR 1
#define SPEX_CHOL_VERSION_MINOR 0
#define SPEX_CHOL_VERSION_SUB   0

// Various other macros are inherited from SPEX_Util.h and may be seen in that file


// Free a pointer and set it to NULL.
#define SPEX_FREE(p)                        \
{                                           \
    SPEX_free (p) ;                         \
    (p) = NULL ;                            \
}                                           \

#ifndef FREE_WORKSPACE
#define FREE_WORKSPACE  
#endif

#define SPEX_MAX(a,b) (((a) > (b)) ? (a) : (b))

#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define SPEX_FLIP(i) (-(i)-2)

#define SPEX_UNFLIP(i) (((i) < 0) ? SPEX_FLIP(i) : (i))

#define SPEX_MARKED(Ap,j) (Ap [j] < 0)

#define SPEX_MARK(Ap,j) { Ap [j] = SPEX_FLIP (Ap [j]) ; }

#define ASSERT(x) assert (x)

// Size of mpz_t, mpq_t and mpfr_t values
#define SIZE_MPZ  sizeof(mpz_t)

#define SIZE_MPQ  sizeof(mpq_t)

#define SIZE_MPFR sizeof(mpfr_t)


// Field access macros for MPZ/MPQ/MPFR struct
// (similar definition in gmp-impl.h and mpfr-impl.h)

#define MPZ_SIZ(x)   ((x)->_mp_size)

#define MPZ_PTR(x)   ((x)->_mp_d)

#define MPZ_ALLOC(x) ((x)->_mp_alloc)

#define MPQ_NUM(x)   mpq_numref(x)

#define MPQ_DEN(x)   mpq_denref(x)

#define MPFR_MANT(x) ((x)->_mpfr_d)

#define MPFR_EXP(x)  ((x)->_mpfr_exp)

#define MPFR_PREC(x) ((x)->_mpfr_prec)

#define MPFR_SIGN(x) ((x)->_mpfr_sign)

#define MPFR_REAL_PTR(x) (&((x)->_mpfr_d[-1])) /*re-define but same result*/

/* Invalid exponent value (to track bugs...) */
#define MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Macros to set the pointer in mpz_t/mpq_t/mpfr_t variable to NULL. It is best
 * practice to call these macros immediately after mpz_t/mpq_t/mpfr_t variable
 * is declared, and before the mp*_init function is called. It would help to
 * prevent error when SPEX_MP*_CLEAR is called before the variable is
 * successfully initialized.
 */

#define SPEX_MPZ_SET_NULL(x)                \
    MPZ_PTR(x) = NULL;                      \
    MPZ_SIZ(x) = 0;                         \
    MPZ_ALLOC(x) = 0;                       \

#define SPEX_MPQ_SET_NULL(x)                \
    MPZ_PTR(MPQ_NUM(x)) = NULL;             \
    MPZ_SIZ(MPQ_NUM(x)) = 0;                \
    MPZ_ALLOC(MPQ_NUM(x)) = 0;              \
    MPZ_PTR(MPQ_DEN(x)) = NULL;             \
    MPZ_SIZ(MPQ_DEN(x)) = 0;                \
    MPZ_ALLOC(MPQ_DEN(x)) = 0;              \

#define SPEX_MPFR_SET_NULL(x)               \
    MPFR_MANT(x) = NULL;                    \
    MPFR_PREC(x) = 0;                       \
    MPFR_SIGN(x) = 1;                       \
    MPFR_EXP(x) = MPFR_EXP_INVALID;         \

/* GMP does not give a mechanism to tell a user when an mpz, mpq, or mpfr
 * item has been cleared; thus, if mp*_clear is called on an object that
 * has already been cleared, gmp will crash. It is also not possible to
 * set a mp*_t = NULL. Thus, this mechanism modifies the internal GMP
 * size of entries to avoid crashing in the case that a mp*_t is cleared
 * multiple times.
 */

#define SPEX_MPZ_CLEAR(x)                   \
{                                           \
    if ((x) != NULL && MPZ_PTR(x) != NULL)  \
    {                                       \
        mpz_clear(x);                       \
        SPEX_MPZ_SET_NULL(x);               \
    }                                       \
}                                           \

#define SPEX_MPQ_CLEAR(x)                   \
{                                           \
    SPEX_MPZ_CLEAR(MPQ_NUM(x));             \
    SPEX_MPZ_CLEAR(MPQ_DEN(x));             \
}                                           \

#define SPEX_MPFR_CLEAR(x)                  \
{                                           \
    if ((x) != NULL && MPFR_MANT(x) != NULL)\
    {                                       \
        mpfr_clear(x);                      \
        SPEX_MPFR_SET_NULL(x);              \
    }                                       \
}                                           \

#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}                                       \

#define SPEX_CHECK(method)              \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        return 0 ;                      \
    }                                   \
}                                       \


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


    
    
//------------------------------------------------------------------------------
// Macros to utilize the default if option is NULL
//------------------------------------------------------------------------------

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
    
static inline int compare (const void * a, const void * b)
{
    return ( *(int64_t*)a - *(int64_t*)b );
}

//------------------------------------------------------------------------------
// SPEX_matrix macros
//------------------------------------------------------------------------------

// These macros simplify the access to entries in a SPEX_matrix.
// The type parameter is one of: mpq, mpz, mpfr, int64, or fp64.

// To access the kth entry in a SPEX_matrix using 1D linear addressing,
// in any matrix kind (CSC, triplet, or dense), in any type:
#define SPEX_1D(A,k,type) ((A)->x.type [k])

// To access the (i,j)th entry in a 2D SPEX_matrix, in any type:
#define SPEX_2D(A,i,j,type) SPEX_1D (A, (i)+(j)*((A)->m), type)

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
// SPEX_Chol_analysis is the data structure for symbolic analysis in the SPEX_Cholesky 
// factorizations. It includes row permutation, elimination tree, and column
// pointers.
//------------------------------------------------------------------------------

typedef struct SPEX_Chol_analysis
{
    int64_t* pinv;      // Row permutation
    int64_t* parent;    // Elimination tree for Cholesky
    int64_t* cp;        // Column pointers for Cholesky
    int64_t lnz;        // Number of nonzeros in Cholesky L
} SPEX_Chol_analysis;
    
    
/* Purpose: Free the Sym_chol data structure */
void SPEX_Chol_analysis_free
(
    SPEX_Chol_analysis* S
);
   
/* Purpose: Permute the matrix A and return A2 = PAP */
SPEX_info SPEX_Chol_permute_A
(
    SPEX_matrix **A2_handle,// Output permuted matrix
    SPEX_matrix* A,        // Initial input matrix
    int64_t* pinv,             // Row permutation
    SPEX_LU_analysis* S    // Column permutation
);

/* Purpose: Compute the elimination tree of A */

int64_t* SPEX_Chol_etree 
(
    SPEX_matrix* A // Input matrix (must be SPD)
);

/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. This is more efficient than the SPEX_reach function 
   It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
int64_t SPEX_Chol_ereach 
(
    SPEX_matrix *A,    // Matrix to be analyzed
    int64_t k,          // Node to start at
    int64_t* parent,    // ELimination Tree
    int64_t* s,         // Contains the nonzero pattern in s[top..n-1]
    int64_t* w          // Workspace array
);

/* Purpose: Depth-first search and postorder of a tree rooted at node j */

int64_t SPEX_Chol_tdfs 
(
    int64_t j,      // Root node
    int64_t k,      
    int64_t* head,  // Head of list
    int64_t* next,  // Next node in the list
    int64_t* post,  // Post ordered tree
    int64_t* stack  // Stack of nodes
);

/* Purpose: post order a forest */
int64_t *SPEX_Chol_post 
(
    int64_t* parent,    // Parent[j] is parent of node j in forest
    int64_t n           // Number of nodes in the forest
);


/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
int64_t SPEX_Chol_leaf 
(
    int64_t i, 
    int64_t j, 
    int64_t* first, 
    int64_t* maxfirst, 
    int64_t* prevleaf,
    int64_t* ancestor, 
    int64_t* jleaf
);

/* Purpose: Something*/
int64_t *SPEX_Chol_counts 
(
    SPEX_matrix *A, 
    int64_t *parent, 
    int64_t *post
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 */
SPEX_info SPEX_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SPEX_matrix* L,              // partial L matrix
    SPEX_matrix* A,              // input matrix
    int64_t k,                    // iteration of algorithm
    int64_t* xi,                  // nonzero pattern vector
    int64_t* parent,              // Elimination tree
    int64_t* c,                   // Column pointers
    SPEX_matrix* rhos,              // sequence of pivots
    int64_t* h,                   // history vector
    SPEX_matrix* x                  // solution of system ==> kth column of L and U
);


/* Purpose: This solves the system L'x = b for Cholesky factorization */
SPEX_info SPEX_Chol_ltsolve 
(
    SPEX_matrix *L,    // The lower triangular matrix
    SPEX_matrix *x      // Solution vector
);

SPEX_info SPEX_Chol_Solve               //solves the linear system LD^(-1)L' x = b
(
    // Output
    SPEX_matrix** x_handle,     // rational solution to the system
    // Input
    SPEX_matrix *A,             // Input matrix (permuted)
    SPEX_matrix* A_orig,        // Input matrix (unpermuted)
    SPEX_matrix* b,             // right hand side vector
    SPEX_matrix* rhos,          // sequence of pivots
    SPEX_matrix* L,             // lower triangular matrix
    int64_t* pinv,                  // row permutation
    SPEX_LU_analysis* S,
    SPEX_options* option        // command options
);

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//--------------------------Alternate Left looking----------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

SPEX_info SPEX_Chol_Factor           // performs the Up lookint64_t*g Cholesky factorization
(
    SPEX_matrix* A,             // matrix to be factored
    SPEX_matrix** L_handle,     // lower triangular matrix
    SPEX_Chol_analysis * S,               // stores guess on nnz and column permutation
    SPEX_matrix ** rhos_handle, // sequence of pivots
    bool left,                  // Set true if performing a left-looking factorization
    SPEX_options* option        // command options
);

/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
SPEX_info SPEX_Chol_Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    SPEX_matrix* A,
    SPEX_matrix** L_handle,              // partial L matrix
    int64_t* xi,                  // nonzero pattern vector
    int64_t* parent,              // Elimination tree
    SPEX_Chol_analysis * S,           // stores guess on nnz and column permutation
    int64_t* c                   // Column pointers
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
SPEX_info SPEX_Left_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SPEX_matrix* L,              // partial L matrix
    SPEX_matrix* A,              // input matrix
    int64_t k,                    // iteration of algorithm
    int64_t* xi,                  // nonzero pattern vector
    SPEX_matrix* rhos,              // sequence of pivots
    int64_t* h,                   // history vector
    SPEX_matrix* x,                  // solution of system ==> kth column of L and U
    int64_t* parent,
    int64_t* c
);

SPEX_info SPEX_Chol_forward_sub
(
    SPEX_matrix *L,   // lower triangular matrix
    SPEX_matrix *x,        // right hand side matrix of size n*numRHS
    SPEX_matrix *rhos      // sequence of pivots used in factorization
);

SPEX_info SPEX_Chol_analyze
(
    SPEX_LU_analysis** S_handle, // symbolic analysis (row/column perm. and nnz L,U)
    const SPEX_matrix *A,        // Input matrix
    const SPEX_options *option   // Control parameters, if NULL, use default
);

SPEX_info SPEX_Chol_backslash
(
    // Output
    SPEX_matrix **X_handle,       // Final solution vector
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_MPQ, SPEX_MPFR, or SPEX_FP64
    const SPEX_matrix *A,         // Input matrix
    const SPEX_matrix *b,         // Right hand side vector(s)
    const SPEX_options* option    // Command options
);




#endif
