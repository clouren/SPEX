//------------------------------------------------------------------------------
// SPEX_Util/Demo/demo.c: Demo of functions in SPEX_Util
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------



// TODO: Do we need a demo in SPEX_Util??
//

#include "../Include/SPEX_Util.h"

#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        printf ("Error: %d line %d file %s\n", ok, __LINE__, __FILE__) ; \
        return 0 ;                      \
    }                                   \
}



// usage:
// example > out
    
int main (void)
{
    
    //--------------------------------------------------------------------------
    // Prior to using SPEX Left LU, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    SPEX_info ok;
    int64_t n = 50, nz = 2500, num=0;
    SPEX_matrix *A = NULL ;                     // input matrix
    SPEX_matrix *R = NULL ;                     // Random matrix to create A
    SPEX_matrix *Rb = NULL;                     // Random matrix to create b
    SPEX_matrix *b = NULL ;                     // Right hand side vector
    SPEX_matrix *x = NULL ;                     // Solution vectors
    SPEX_options *option = SPEX_create_default_options();
    if (!option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        //FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Generate a random dense 50*50 matrix
    //--------------------------------------------------------------------------

    // R is a n*n triplet matrix whose entries are FP64 Note that the first
    // boolean parameter says that the matrix is not shallow, so that A->i,
    // A->j, and A->x are calloc'd. The second boolean parameter is meaningless
    // for FP64 matrices, but it tells SPEX LU to allocate the values of A->x
    // for the mpz_t, mpq_t, and mpfr_t entries
    SPEX_matrix_allocate(&R, SPEX_TRIPLET, SPEX_FP64, n, n, nz,
        false, true, option);
    
    // Rb is a n*1 dense matrix whose entries are FP64
    SPEX_matrix_allocate(&Rb, SPEX_DENSE, SPEX_FP64, n, 1, n,
        false, true, option);

    // Randomly generate the input
    unsigned int seed = 10;
    srand(seed);
    for (int64_t k = 0; k < n; k++)
    {
        Rb->x.fp64[k] = rand();
        for (int64_t p = 0; p < n; p++)
        {
            R->i[num] = k;
            R->j[num] = p;
            R->x.fp64[num] = rand();
            num+=1;
        }
    }

    R->nz = n*n;

    //--------------------------------------------------------------------------
    // Build A and b
    //--------------------------------------------------------------------------

    // A is a copy of the R matrix. A is a CSC matrix with mpz_t entries
    OK ( SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, R, option));
    // b is a copy of the Rb matrix. b is dense with mpz_t entries. 
    OK ( SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, Rb, option));

    int64_t nnz2 = SPEX_matrix_nnz(A,option);
    printf("\nnnz2 is: %ld\n",nnz2);
    
    printf("\n\nHello\n\n");
    return 0;
}

