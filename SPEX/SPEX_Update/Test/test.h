//------------------------------------------------------------------------------
// SPEX_Update/Test/test.h: #include file for the test programs
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2023, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------
#include <ctype.h>
#include <glpk.h>
#include "SPEX.h"

#ifndef FREE_WORKSPACE
#define FREE_WORKSPACE
#endif

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

#define OK(method)                               \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_OK)                         \
    {                                            \
        FREE_WORKSPACE;                          \
        SPEX_PRINT_INFO(info);                   \
        return 0 ;/*continue; */                 \
    }                                            \
}
// These flags are used for code development.  Uncomment them as needed.

// to turn on debugging, uncomment this line:
// #undef NDEBUG

#undef ASSERT

#ifndef NDEBUG

    // debugging enabled
    #ifdef MATLAB_MEX_FILE
    #define ASSERT(x) \
    {                                                                       \
        if (!(x))                                                           \
        {                                                                   \
            mexErrMsgTxt ("failure: " __FILE__ " line: "                    \
                LAGRAPH_XSTR(__LINE__)) ;                                   \
        }                                                                   \
    }
    #else
    #include <assert.h>
    #define ASSERT(x) assert (x) ;
    #endif

#else

    // debugging disabled
    #define ASSERT(x)

#endif

//------------------------------------------------------------------------------
// Matrix Market format
//------------------------------------------------------------------------------

// %%MatrixMarket matrix <fmt> <type> <storage> uses the following enums:

typedef enum
{
    MM_coordinate,
    MM_array,
}
MM_fmt_enum ;

typedef enum
{
    MM_real,
    MM_integer,
    MM_complex,
    MM_pattern
}
MM_type_enum ;

typedef enum
{
    MM_general,
    MM_symmetric,
    MM_skew_symmetric,
    MM_hermitian
}
MM_storage_enum ;

// maximum length of each line in the Matrix Market file format

// The MatrixMarket format specificies a maximum line length of 1024.
// This is currently sufficient for GraphBLAS but will need to be relaxed
// if this function is extended to handle arbitrary user-defined types.
#define MMLEN 1024
#define MAXLINE MMLEN+6

SPEX_info SPEX_mmread
(
    SPEX_matrix *A_handle,// handle of matrix to create
    FILE *f,             // file to read from, already open
    SPEX_options option
);

// read the remaining component (b, lb, ub and c) to finish constructing LP
SPEX_info SPEX_construct_LP
(
    glp_prob *LP,
    SPEX_matrix *A_handle,
    SPEX_matrix *b_handle,
    SPEX_matrix *c_handle,
    double *z0_handle,
    char *file_name,
    SPEX_options option
);


SPEX_info SPEX_A_plus_vvT
(
    SPEX_matrix A0,
    const SPEX_matrix M,
    const int64_t j
);

SPEX_info SPEX_matrix_equal
(
    bool *Isequal,
    const SPEX_matrix L1,
    const SPEX_matrix L_update,
    const int64_t *P_update
);

SPEX_info MY_update_verify
(
    bool *Is_correct,     // if the factorization is correct
    SPEX_factorization F,// LU factorization of A
    const SPEX_matrix A,     // Input matrix
    const SPEX_options option// command options
);
