//------------------------------------------------------------------------------
// SPEX_WAMF/Demo/example.c: example main program for SPEX_WAMF
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------


/* This example shows how to use SPEX WAMF with a given input matrix 
 * The input is read from a file */

// usage:
// example > out
// out is file for output calculated result

#define FREE_WORKSPACE              \
    SPEX_LU_analysis_free(&S, option);\
    SPEX_matrix_free(&A, option);   \
    SPEX_FREE(option);              \
    SPEX_finalize();

#include "demos.h"   

#define DEMO_OK(method)                 \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}

    
int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX WAMF, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Get matrix and right hand side file names
    //--------------------------------------------------------------------------
    char *mat_name, *rhs_name;
    mat_name = "../ExampleMats/10teams_mat.txt";
    rhs_name = "../ExampleMats/10teams_mat.soln.txt";
    if (argc > 2)
    {
        mat_name = argv[1];
        rhs_name = argv[2];
    }

    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_info ok;
    SPEX_matrix *A = NULL ;                     // input matrix
    SPEX_LU_analysis *S = NULL ;                // Column permutation (COLAMD)
    SPEX_LU_analysis *S2 = NULL;                // Column permutation (AMD)
    
    SPEX_matrix *B = NULL;                      // Bit-matrix copy of A
    
    SPEX_options *option = NULL;
    SPEX_create_default_options(&option);
    if (option == NULL)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Allocate memory, read in A
    //--------------------------------------------------------------------------

    // Read in A. The output of this demo function is A in CSC format with
    // mpz_t entries.
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    DEMO_OK(SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);
    
        
    //--------------------------------------------------------------------------
    // Create copy of A
    //--------------------------------------------------------------------------
    B = SPEX_WAMF_get_bitmat(A);
    
    //--------------------------------------------------------------------------
    // Perform Orderings of A
    //--------------------------------------------------------------------------
    
    S = (SPEX_LU_analysis*) SPEX_malloc(sizeof(SPEX_LU_analysis));

    // Allocate memory for column permutation
    S->q = (int64_t*) SPEX_malloc((A->n+1) * sizeof(int64_t));
    
    S2 = (SPEX_LU_analysis*) SPEX_malloc(sizeof(SPEX_LU_analysis));

    // Allocate memory for column permutation
    S2->q = (int64_t*) SPEX_malloc((A->n+1) * sizeof(int64_t));
    
    clock_t start_col = clock();
    option->order = SPEX_AMD;  // AMD
    double Control [AMD_CONTROL];           // Declare AMD control
    amd_l_defaults (Control) ;              // Set AMD defaults
    double Info [AMD_INFO];
    // Perform AMD
    amd_l_order(A->n, (SuiteSparse_long *) A->p, (SuiteSparse_long *) A->i,
            (SuiteSparse_long *) S->q, Control, Info) ;
    S->lnz = Info[AMD_LNZ];        // estimate for unz and lnz    
    clock_t end_col = clock();
    
    clock_t start_col2 = clock();
    option->order = SPEX_COLAMD;  // COLAMD
    // Declared as per COLAMD documentation
    int64_t Alen = 2*A->nzmax + 6 *(A->n+1) + 6*(A->n+1) + A->n;
    int64_t* A2 = (int64_t*) SPEX_malloc(Alen* sizeof(int64_t));
    if (!A2)
    {
        // out of memory
        SPEX_LU_analysis_free (&S, option) ;
        return (SPEX_OUT_OF_MEMORY) ;
    }
    // Initialize S->q as per COLAMD documentation
    for (int64_t i = 0; i < A->n+1; i++)
    {
        S2->q[i] = A->p[i];
    }
    // Initialize A2 per COLAMD documentation
    for (int64_t i = 0; i < A->nzmax; i++)
    {
        A2[i] = A->i[i];
    }
    int64_t stats [COLAMD_STATS];
    colamd_l (A->n, A->n, Alen, (SuiteSparse_long *) A2,
            (SuiteSparse_long *) S2->q, (double *) NULL,
            (SuiteSparse_long *) stats) ;
    // estimate for lnz and unz
    S2->lnz = 10*A->nzmax;    
    clock_t end_col2 = clock();
    
    clock_t start_col3 = clock();
    
    // How A is analyzed. If 0 we analyze A directly. If 1 we analyze A'+A. If 2 we analyze A'*A
    int64_t order = 1;
    // How the weights are computed. If 0 weighted sum, if 1 sum, if 2 mean, if 3 min,
    // if 4 max, if 5 diagonal
    int64_t option2 = 0;
    // Convex combination parameter. Larger it is the more sparsity is favored
    double alpha = 0.8;
    // True if eliminate dense columns prior to factorization
    bool aggressive = true;
    
    int64_t* q1 = SPEX_WAMF_wamf(order, B, option2, alpha, aggressive);
    if (q1 == NULL)
        printf("\nHELLO!");
    clock_t end_col3 = clock();
    
    
    clock_t start_col4 = clock();
    // How A is analyzed. If 0 we analyze A directly. If 1 we analyze A'+A. If 2 we analyze A'*A
    order = 2;
    int64_t* q2 = SPEX_WAMF_wamf(order, B, option2, alpha, aggressive);
    if (q2 == NULL)
        printf("\nHEY HEY HEY!");
    clock_t end_col4 = clock();
    
    
    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    
    int64_t k;
    printf("\nAMD is:\n");
    for (k = 0; k < A->n; k++)
        printf("%ld ",S->q[k]);
    printf("\nCOLAMD is:\n");
    for (k = 0; k < A->n; k++)
        printf("%ld ",S2->q[k]);
    printf("\nWAMF A'+A is:\n");
    for (k = 0; k < A->n; k++)
        printf("%ld ",q1[k]);
    printf("\nWAMF A'*A is:\n");
    for (k = 0; k < A->n; k++)
        printf("%ld ",q2[k]);
    
    double t_col = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_col2 = (double) (end_col2-start_col2)/CLOCKS_PER_SEC;
    double t_col3 = (double) (end_col3-start_col3)/CLOCKS_PER_SEC;
    double t_col4 = (double) (end_col4-start_col4)/CLOCKS_PER_SEC;
    
    printf("\nAMD time: \t\t\t\t%lf", t_col);
    printf("\nCOLAMD time: \t\t\t\t%lf", t_col2);
    printf("\nWAMF A'+A time: \t\t\t%lf", t_col3);
    printf("\nWAMF A'*A time: \t\t\t%lf", t_col4);
    
    
    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;

    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}

