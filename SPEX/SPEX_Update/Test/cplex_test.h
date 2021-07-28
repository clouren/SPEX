//------------------------------------------------------------------------------
// SPEX_Update/Test/test.h: #include file for the test programs
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------
#include <ctype.h>
#include <glpk.h>
#include <ilcplex/cplex.h>
#include "SPEX_Update.h"
#include "SPEX_Left_LU.h"

#ifndef FREE_WORKSPACE
#define FREE_WORKSPACE
#endif

#define GOTCHA \
    printf ("%s, line %d\n", __FILE__, __LINE__);

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
    SPEX_matrix **A_handle,// handle of matrix to create
    FILE *f,             // file to read from, already open
    SPEX_options *option
);

// read the remaining component (b, lb, ub and c) to finish constructing LP
SPEX_info SPEX_construct_LP
(
    glp_prob *LP,
    SPEX_matrix **A_handle,
    SPEX_matrix **b_handle,
    SPEX_matrix **c_handle,
    double *z0_handle,
    char *file_name,
    SPEX_options *option
);

SPEX_info SPEX_get_CPLEX_LP
(
    CPXENVptr env,
    CPXLPptr *LP,
    int *status,
    SPEX_matrix **A_handle,
    SPEX_matrix **b_handle,
    SPEX_matrix **c_handle,
    char *file_name,
    SPEX_options *option
);
