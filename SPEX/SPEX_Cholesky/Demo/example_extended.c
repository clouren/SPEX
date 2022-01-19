//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol: Demo main program for SPEX_Chol
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Include the Integer-preserving Cholesky routines */

#define FREE_WORKSPACE                  \
{                                       \
    SPEX_matrix_free(&A,NULL);          \
    SPEX_matrix_free(&PAP,NULL);         \
    SPEX_matrix_free(&b,NULL);          \
    SPEX_matrix_free(&x,NULL);          \
    SPEX_symbolic_analysis_free (&S, option); \
    SPEX_factorization_free(&F, option);     \
    SPEX_FREE(option);                  \
    SPEX_finalize();                    \
}                                       
    //SPEX_symbolic_analysis_free (&S, option); \ //TODO double check because not everything in S is initialized because we don't need everything that gets freed
    //TOASK SPEX_factorization_free(&F, option) does not seem to be freeing rhos and L ???
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
   
    SPEX_symbolic_analysis *S = NULL;
    SPEX_factorization *F = NULL ;
    SPEX_matrix* A = NULL;
    SPEX_matrix* b = NULL;
    SPEX_matrix* PAP = NULL;
    SPEX_matrix* x = NULL;
    
    // Default options. May be changed in SLIP_LU_config.h
    SPEX_options *option = NULL;
    DEMO_OK(SPEX_create_default_options(&option));
    
    char* mat_name = "../ExampleMats/2.mat.txt";// Set demo matrix and RHS name
    char* rhs_name = "../ExampleMats/2.mat.soln.txt";
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
    
    DEMO_OK(SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);
    n = A->n;
    // For this code, we utilize a vector of all ones as the RHS vector    
    SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true, option);
    // Create RHS
    for (int64_t k = 0; k < n; k++)
        DEMO_OK(SPEX_mpz_set_ui(b->x.mpz[k],1));
    
    //--------------------------------------------------------------------------
    // Perform Ordering of A
    //--------------------------------------------------------------------------
    clock_t start_col = clock();
        
    // Symmetric ordering of A. Uncomment the desired one, AMD is recommended
    //option->order = SPEX_NO_ORDERING;  // No ordering
    option->order = SPEX_AMD;  // AMD
    //option->order = SPEX_COLAMD; // COLAMD

    DEMO_OK(SPEX_Chol_preorder(&S, A, option));    
    clock_t end_col = clock();
    
    
    //--------------------------------------------------------------------------
    // Determine if A is indeed symmetric. If so, we try Cholesky
    // uncomment the one desired.
    // --------------------------------------------------------------------------
    
    clock_t start_sym = clock();
    bool test;
    //DEMO_OK(SPEX_determine_symmetry(A, 0, option));    // Determine symmetry just with nonzero pattern
    DEMO_OK( SPEX_determine_symmetry(A, 1, option));    // Determine symmetry with nonzero pattern and values
   
    clock_t end_sym = clock();
    //--------------------------------------------------------------------------
    // Permute matrix A, that is set PAP = PAP'
    //--------------------------------------------------------------------------
    /*
    pinv = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv[index] = k;
    }
    */
    
    DEMO_OK( SPEX_Chol_permute_A(&PAP, A, S));
    option->print_level = 3;
    option->check = true;
    
  //SPEX_matrix_check(PAP,option);

    //--------------------------------------------------------------------------
    // SPEX Chol Factorization
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();
    
    
    SPEX_matrix* L2 = NULL;
    SPEX_matrix* rhos2 = NULL;
    DEMO_OK( SPEX_Chol_Factor(&F, S,PAP, option));
    
//    SPEX_matrix_check(L, option);
//     
     F->L->m = n;
//     
     clock_t end_factor = clock();
     
     //printf("\nhere\n");

    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();
    option->check = true;

    DEMO_OK( SPEX_Chol_Solve(&x, F, b, option));
    
    clock_t end_solve = clock();

    if (option->check)
    {
        DEMO_OK(SPEX_check_solution(A, x, b, option));
    }
    
    
    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------
    
    double t_col = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_sym = (double) (end_sym-start_sym)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;

    printf("\nNumber of L nonzeros: \t\t\t%ld",
        (F->L->p[F->L->n]) );
    printf("\nSymmetry Check time: \t\t\t%lf", t_sym);
    printf("\nSymbolic analysis time: \t\t%lf", t_col);
    printf("\nIP Chol Factorization time: \t\t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);
    

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
                 
}
    
