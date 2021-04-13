//------------------------------------------------------------------------------
// SPEX/Tests/SPEX_WAMF_Unsymmetric: Test SPEX Left with WAMF
//------------------------------------------------------------------------------

// SPEX: (c) 2021

//------------------------------------------------------------------------------

/* Include the Integer-preserving Cholesky routines */

#define FREE_WORKSPACE                  \
{                                       \
    SPEX_matrix_free(&A,NULL);          \
    SPEX_matrix_free(&A2,NULL);         \
    SPEX_matrix_free(&b,NULL);          \
    SPEX_matrix_free(&x,NULL);          \
    SPEX_FREE(q1)                       \
    SPEX_Chol_analysis_free(&S);        \
    SPEX_FREE(option);                  \
    SPEX_finalize();                    \
}                                       \

#define DEMO_OK(method)                 \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        SPEX_Chol_determine_error(ok);  \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}


# include "demos.h"

/* The purpose of this function is to compare SPEX Left LU
 * with COLAMD vs SPEX Left LU with WAMF on BasisLIB matrices. 
 * A similar file exists for SJ Database matrices and 
 * for Cholesky factorization
 */


int main( int argc, char* argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX-Chol, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function. 
    //--------------------------------------------------------------------------
    SPEX_initialize();
    
    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int64_t n = 0, check, ok, j, index, k, nz = 0;
    
    SPEX_matrix *A = NULL;
    SPEX_matrix *A2 = NULL;
    SPEX_matrix *b = NULL;
    SPEX_matrix* x = NULL;
    SPEX_Chol_analysis* S = NULL;
    int64_t* q1 = NULL;
    
    // Default options. May be changed in SLIP_LU_config.h
    SPEX_options *option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));
    
    char* mat_name = "10teams_mat.txt"; // Set demo matrix and RHS name
    char* rhs_name = "10teams_v.txt";
    int64_t rat = 1;
    
    // Process the command line
    DEMO_OK(SPEX_Chol_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));
    
    //--------------------------------------------------------------------------
    // Allocate memory    
    //--------------------------------------------------------------------------
    
    // Read in A
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    
    DEMO_OK(SPEX_tripread(&A, mat_file, option));
    fclose(mat_file);
    n = A->n;
    
    FILE* rhs_file = fopen(rhs_name,"r");
    if( rhs_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    
    DEMO_OK(SPEX_read_dense(&b, rhs_file, option));
    fclose(rhs_file);

    //--------------------------------------------------------------------------
    // Solve Ax = b with COLAMD
    //--------------------------------------------------------------------------
    
    clock_t start_s = clock();
   
    // SPEX Left LU has an optional check, to enable it, one can set the following
    // parameter to be true.
    option->order = SPEX_COLAMD;
    //option->check = true;
    // Solve the system and give double solution
    OK(SPEX_Left_LU_backslash( &x, SPEX_FP64, A, b, option));
         
    clock_t end_s = clock();

    double t_s = (double) (end_s - start_s) / CLOCKS_PER_SEC;

    printf("\nSPEX Left LU COLAMD Factor & Solve time: %lf\n", t_s);
    
    //--------------------------------------------------------------------------
    // Do WAMF on A
    //--------------------------------------------------------------------------
    
    clock_t start_w = clock();
    option->order = SPEX_NO_ORDERING;
    
    // Create copy of A
    SPEX_matrix*B = SPEX_WAMF_get_bitmat(A);
    
    // Perform WAMF
    
    // How A is analyzed. If 0 we analyze A directly. If 1 we analyze A'+A. If 2 we analyze A'*A
    int64_t W_order = 2;
    // How the weights are computed. If 0 weighted sum, if 1 sum, if 2 mean, if 3 min,
    // if 4 max, if 5 diagonal
    int64_t W_option2 = 0;
    // Convex combination parameter. Larger it is the more sparsity is favored
    double W_alpha = 0.8;
    // True if eliminate dense columns prior to factorization
    bool W_aggressive = true;
    
    q1 = SPEX_WAMF_wamf(W_order, B, W_option2, W_alpha, W_aggressive);
    
    
    S = (SPEX_Chol_analysis*) SPEX_malloc(sizeof(SPEX_Chol_analysis));
    S->q = (int64_t*) SPEX_malloc((n+1) * sizeof(int64_t));
    //pinv = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    for (k = 0; k < n; k++)
    {
        S->q[k] = q1[k];
    }
    
    
    DEMO_OK( SPEX_Chol_permute_A(&A2, A, S));
    
    OK(SPEX_Left_LU_backslash( &x, SPEX_FP64, A2, b, option));
         
    clock_t end_w = clock();

    double t_w = (double) (end_w - start_w) / CLOCKS_PER_SEC;

    printf("\nSPEX Left LU WAMF Factor & Solve time: %lf\n", t_w);
    
    
/*   
    
    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------
    
    printf("\nNumber of L nonzeros: \t\t\t%ld",
        (L->p[L->n]) );
    printf("\nSymmetry Check time: \t\t\t%lf", t_sym);
    printf("\nSymbolic analysis time: \t\t%lf", t_col);
    printf("\nIP Chol Factorization time: \t\t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);
*/
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
}
    
