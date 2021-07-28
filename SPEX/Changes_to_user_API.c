// data structure for factorization
typedef enum
{
    SPEX_LU_FACTORIZATION = 0,            // LU factorization
    SPEX_LU_ANALYSIS = 1,                 // LU analysis
    SPEX_CHOLESKY_FACTORIZATION = 2,      // Cholesky factorization
    SPEX_CHOLESKY_ANALYSIS = 3,           // Cholesky analysis
    SPEX_QR_FACTORIZATION = 4,            // QR factorization
    SPEX_QR_ANALYSIS = 5                  // QR factorization
}SPEX_factorization_kind ;

typedef struct
{
    SPEX_factorization_type kind;
    SPEX_matrix *L; // check if dynamic from L->kind
    SPEX_matrix *U;
    SPEX_matrix *Q;
    SPEX_matrix *R;
    SPEX_matrix *rhos;
    int64_t *P_perm;// should not free by Left_LU_factorize
    int64_t *Pinv_perm;
    int64_t *Q_perm;// from SPEX_*_analysis *S
    int64_t *Qinv_perm;
    mpq_t scale_for_A;
    int64_t lnz;// consider combine with SPEX_*_analysis
    int64_t unz;
    int64_t* parent;    // Elimination tree of A for Cholesky
    int64_t* cp;        // Column pointers of L
} SPEX_factorization;

// for reference
#if 0
typedef struct SPEX_Chol_analysis
{
    int64_t* pinv;      // Row permutation
    int64_t* q ;        // Column permutation, representing
                        // the permutation matrix P. The matrix P*A*P' is factorized.
                        // If the kth column of L, and P*A*P' is column j of the
                        // unpermuted matrix A, then j = S->q [k].
    int64_t* parent;    // Elimination tree of A for Cholesky
    int64_t* cp;        // Column pointers of L
    int64_t lnz;        // Number of nonzeros in Cholesky L (might be estimate)
                        // at initialization if default column ordering (AMD) is used
                        // it will be exact otherwise it will be an estimate
                        // after elimination tree is computed it will be exact
} SPEX_Chol_analysis;
#endif

//------------------------------------------------------------------------------
// functions in SPEX_Util
//------------------------------------------------------------------------------
// TODO move SPEX_LU_analyze etc. to SPEX_Left_LU/Source

SPEX_info SPEX_analyze
{
   SPEX_factorization_kind kind,
   SPEX_factorization *S,
   const SPEX_options* option // command options
}

//------------------------------------------------------------------------------
// functions in SPEX_Left_LU
//------------------------------------------------------------------------------
// old
SPEX_info SPEX_Left_LU_factorize
(
    // output:
    SPEX_matrix **L_handle,    // lower triangular matrix
    SPEX_matrix **U_handle,    // upper triangular matrix
    SPEX_matrix **rhos_handle, // sequence of pivots
    int64_t **pinv_handle,     // inverse row permutation
    // input:
    const SPEX_matrix *A,      // matrix to be factored
    const SPEX_LU_analysis *S, // column permutation and estimates
                               // of nnz in L and U 
    const SPEX_options* option // command options
);

// new
SPEX_info SPEX_factorize
(
    // output:
    SPEX_factorization **F,    // LU factorization of A
    // input:
    const SPEX_matrix *A,      // matrix to be factored
    const SPEX_factorization *S,// S->kind = SPEX_LU_analysis
    const SPEX_options* option // command options
);

/*
use case: if A1 and A1 are two different matrices with same patter, we can do
          SPEX_LU_analyze(&S, A1, opt);
          SPEX_Left_LU_factorize(&F1, A1, S, opt);
          SPEX_Left_LU_factorize(&F2, A2, S, opt);
          ...
*/

//------------------------------------------------------------------------------
// old
SPEX_info SPEX_Left_LU_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_matrix **x_handle,  // rational solution to the system
    // input:
    const SPEX_matrix *b,   // right hand side vector
    const SPEX_matrix *A,   // Input matrix
    const SPEX_matrix *L,   // lower triangular matrix
    const SPEX_matrix *U,   // upper triangular matrix
    const SPEX_matrix *rhos,// sequence of pivots
    const SPEX_LU_analysis *S,// symbolic analysis struct
    const int64_t *pinv,    // inverse row permutation
    const SPEX_options* option // Command options
);

//new
SPEX_info SPEX_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_matrix **x_handle,  // rational solution to the system
    // input:
    const SPEX_matrix *b,   // right hand side vector
    // const SPEX_matrix *A,   // Input matrix
    // TODO A is now only needed for SPEX_check_solution, which could be
    // modified and isolated from SPEX_Left_LU_solve to remove the parameter
    // "A".  Refer to SPEX_update_verify
    const SPEX_factorization *F,// upper triangular matrix
    const SPEX_options* option // Command options
);
// TODO enable to let user call with SPEX_Left_LU_solve(&b, NULL, F, option)?
// Refer to SPEX_Update_Solve

//------------------------------------------------------------------------------
// functions in SPEX_Update
//------------------------------------------------------------------------------
// old
SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_matrix *L,   // n-by-n dynamic_CSC matrix that gives the Cholesky
                      // factorization
    SPEX_matrix *rhos,// n-by-1 dense matrix that gives the array of pivots
    const int64_t *P, // row permutation
    const int64_t *P_inv,// inverse of row permutation
    SPEX_vector *w,   // a n-by-1 vector in sparse compressed column form that
                      // modifies the original matrix A, the resulting A is
                      // A+sigma*w*w^T. In output, w is updated as the solution
                      // to L*D^(-1)*w_out = w
    const int64_t sigma,// a scalar that determines whether this is an update
                      // or downdate
    const SPEX_options *option
);

//new
SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_factorization *F,// dynamic cholesky factorization
    SPEX_vector *w,   // a n-by-1 vector in sparse compressed column form that
                      // modifies the original matrix A, the resulting A is
                      // A+sigma*w*w^T. In output, w is updated as the solution
                      // to L*D^(-1)*w_out = w
    const int64_t sigma,// a scalar that determines whether this is an update
                      // or downdate
    const SPEX_options *option
);


//------------------------------------------------------------------------------
// old    
SPEX_info SPEX_Update_LU_ColRep 
(   
    SPEX_matrix *A,         // n-by-n original matrix in dynamic_CSC form 
    SPEX_matrix *L,         // n-by-n lower triangular factorization of A in 
                            // dynamic_CSC form.
    SPEX_matrix *UT,        // The transpose of U in dynamic_CSC form, where U 
                            // is the n-by-n upper triangular factorization of 
                            // A.
    SPEX_matrix *rhos,      // n-by-1 dense matrix for the array of pivots 
    int64_t *P,             // Row permutation 
    int64_t *P_inv,         // Inverse of row permutation 
    int64_t *Q,             // Column permutation 
    int64_t *Q_inv,         // Inverse of column permutation 
    SPEX_vector **vk,       // Pointer to the inserted column in the compressed 
                            // column form, the rows of vk are in the same orde 
                            // as A. This vector will be swapped with A->v[k] 
                            // in the output upon return, regardless of failure.
    int64_t k,              // The column index that vk will be inserted 
    const SPEX_options *option// Command parameters 
);


//new    
SPEX_info SPEX_Update_LU_ColRep 
(   
    SPEX_matrix *A,         // n-by-n original matrix in dynamic_CSC form 
    SPEX_factorization *F,  // dynamic LU factorization
    SPEX_vector **vk,       // Pointer to the inserted column in the compressed 
                            // column form, the rows of vk are in the same order
                            // as A. This vector will be swapped with A->v[k] 
                            // in the output upon return, regardless of failure 
    int64_t k,              // The column index that vk will be inserted 
    const SPEX_options *option// Command parameters 
);


//------------------------------------------------------------------------------
// old
SPEX_info SPEX_Update_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a n*m dense matrix contains the solution to
                            // the system. If users wish to overwrite the
                            // solution to the right-hand-side matrix b, this
                            // can be provided as &b. Otherwise, new space will
                            // be allocated for x_handle
    // input:
    SPEX_matrix *b,         // a n*m dense matrix contains the right hand
                            // side vector
    const SPEX_matrix *L,   // a n*n dynamic_CSC matrix that gives the lower
                            // triangular matrix
    const SPEX_matrix *UT,  // a n*n dynamic_CSC matrix that gives the transpose
                            // of the upper triangular matrix
    const mpq_t A_scale,    // scale of the input matrix
    int64_t *h,             // history vector
    const SPEX_matrix *rhos,// a n*1 dense matrix that gives the array of pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv,   // inverse of column permutation
    const SPEX_options* option // Command options
);


//new
SPEX_info SPEX_Update_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a n*m dense matrix contains the solution to
                            // the system. If users wish to overwrite the
                            // solution to the right-hand-side matrix b, this
                            // can be provided as &b. Otherwise, new space will
                            // be allocated for x_handle
    // input:
    SPEX_matrix *b,         // a n*m dense matrix contains the right hand
                            // side vector
    const SPEX_factorization *F,// a dynamic factorization (either LU or chol)
    int64_t *h,             // history vector
    const SPEX_options* option // Command options
);
