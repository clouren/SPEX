//------------------------------------------------------------------------------
// IP-Chol/MATLAB/IP_Chol_mex_soln: Use IP-Chol within MATLAB
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "./Source/IP-Chol_mex.h"

/* Purpose: .c files defining the IP-Chol matlab interfacee
 * This function defines: x = IP_Chol(A, b, option)
 */

//TODO Mimic SPEX_backslash
//TODO Fix MATLAB interface

void mexFunction
(
    int32_t nargout,
    mxArray *pargout [ ],
    int32_t nargin,
    const mxArray *pargin [ ]
)
{
    //--------------------------------------------------------------------------
    // Initialize SPEX LU library environment
    //--------------------------------------------------------------------------
    SPEX_initialize_expert (mxMalloc, mxCalloc, mxRealloc, mxFree) ;
    SPEX_info status;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    IP_check_input(pargin, nargin);
    if (nargout > 1 || nargout <= 0 || nargin != 3)
    {
        mexErrMsgTxt("Usage: x = IP_Chol(A,b,option)");
    }
    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------
    SPEX_LU_analysis* S = NULL;
    SPEX_matrix *A = NULL;
    SPEX_matrix *L = NULL;
    SPEX_matrix *b = NULL;
    SPEX_matrix *rhos = NULL;
    int64_t* pinv = NULL;
    int64_t* pinv2 = NULL;
    SPEX_matrix* A2 = NULL;
    Sym_chol* S2 = NULL;
    SPEX_matrix* x = NULL;
    SPEX_options *option = SPEX_create_default_options();
    //if (!A || !L || !b || !option)
   // {
    //    IP_mex_error (SPEX_OUT_OF_MEMORY);
   // }

    //--------------------------------------------------------------------------
    // Declare variables and process input
    //--------------------------------------------------------------------------
    // Read in options
    IP_get_matlab_options(option, pargin[2]);

    // Read in A and b
    IP_mex_get_A_and_b(&A, &b, pargin, nargin, option);

    //--------------------------------------------------------------------------
    // Perform Ordering of A
    //--------------------------------------------------------------------------        
    IP_MEX_OK(SPEX_LU_analyze(&S, A, option));    
    
    //--------------------------------------------------------------------------
    // Permute matrix A, that is set A2 = PAP'
    //--------------------------------------------------------------------------
    int64_t n = A->n;
    pinv2 = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    for (int64_t k = 0; k < n; k++)
    {
        int64_t index = S->q[k];
        pinv2[index] = k;
    }
    
    IP_MEX_OK(IP_Chol_permute_A(&A2, A, pinv2, S));
    
    
    
    //--------------------------------------------------------------------------
    // SPEX Chol Factorization
    //--------------------------------------------------------------------------
    
    S2 = (Sym_chol*) SPEX_malloc(1* sizeof(Sym_chol));
    
    IP_MEX_OK(IP_Chol_Factor( A2, &L, S2, &rhos, false, option));
    
    
    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    
    IP_MEX_OK(IP_Solve( &x, A2, A, b, rhos, L, pinv2, S, option));
    
    SPEX_matrix *x2 = NULL;
    IP_MEX_OK(SPEX_matrix_copy(&x2, SPEX_DENSE, SPEX_FP64, x, option));
    
    mxArray* Xmatlab = mxCreateDoubleMatrix ((mwSize) x2->m, (mwSize) x2->n,
        mxREAL);
    double* x_out = mxGetPr(Xmatlab);

    for (int k = 0; k < x->n*x->m; k++)
    {
        x_out[k] = x2->x.fp64[k];
    }
    
    pargout[0] =  Xmatlab;
    SPEX_matrix_free(&b, option);
    SPEX_matrix_free(&A, option);
    SPEX_matrix_free(&A2, option);
    SPEX_FREE(option);
    SPEX_finalize();
}
