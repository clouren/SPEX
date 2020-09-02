//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Solve: Solve the SPD linear system after factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE                      \
    SPEX_matrix_free(&b2, option);          \
    SPEX_matrix_free(&x, option);           \
    SPEX_MPQ_CLEAR(scale);                  \
    SPEX_MPQ_CLEAR(det2);                   \

#include "SPEX_Chol.h"
    
/* Purpose: This function solves the linear system LD^(-1)L' x = b.*
 *
 * Input arguments:
 * 
 * x_handle:        A handle to the solution matrix. On input this is NULL,
 *                  on output x_handle contains a pointer to the solution vector(s)
 * 
 * A:               Permuted version of the input matrix
 * 
 * A_orig:          Nonpermuted input matrix
 * 
 * b:               Right hand side vector(s)
 * 
 * rhos:            Sequence of pivots encountered in factorization
 * 
 * L:               Lower triangular matrix
 * 
 * pinv:            Inverse row permutation
 * 
 * S:               Column permutation struct
 * 
 * option:          Command options
 * 
 */
SPEX_info SPEX_Chol_Solve       // solves the linear system LDL' x = b
(
    // Output
    SPEX_matrix** x_handle,     // rational solution to the system
    // Input
    SPEX_matrix* A,             // Input matrix (permuted)
    SPEX_matrix* A_orig,        // Nonpermuted input matrix
    SPEX_matrix* b,             // right hand side vector
    SPEX_matrix* rhos,          // sequence of pivots
    SPEX_matrix* L,             // lower triangular matrix
    int64_t* pinv,              // row permutation
    SPEX_LU_analysis *S,        // Column permutation
    SPEX_options* option        // command options
)
{
    SPEX_info ok;
    // Check the inputs
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    SPEX_REQUIRE(A_orig, SPEX_CSC, SPEX_MPZ);
    SPEX_REQUIRE(b, SPEX_DENSE, SPEX_MPZ);
    SPEX_REQUIRE(rhos, SPEX_DENSE, SPEX_MPZ);
    SPEX_REQUIRE(L, SPEX_CSC, SPEX_MPZ);
    
    if (!x_handle || !pinv || !S) return SPEX_INCORRECT_INPUT;
    
    int64_t i, j, k, n = L->n, nz;
    
    mpq_t scale, det2 ;
    SPEX_MPQ_SET_NULL (scale) ;
    SPEX_MPQ_SET_NULL (det2) ;

    SPEX_matrix *x = NULL;   // unpermuted solution
    SPEX_matrix *x2 = NULL;  // permuted final solution
    SPEX_matrix *b2 = NULL;  // permuted b
    
    // Permute b and place it in b2
    SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPZ, b->m, b->n, b->m*b->n, false, true, option);
    for (i = 0; i < b->m; i++)
    {
        for (j = 0; j < b->n; j++)   
        {
            OK (SPEX_mpz_set( SPEX_2D(b2, pinv[i], j, mpz),
                              SPEX_2D(b, i, j, mpz)));
        }
    }
    
    // b2 = L \ b2    
    OK(SPEX_Chol_forward_sub(L, b2, rhos));
    
    // b2 = b2 * det 
    nz = SPEX_matrix_nnz(b, NULL);
    for (i = 0; i < nz; i++)
    {
        OK ( SPEX_mpz_mul( b2->x.mpz[i], b2->x.mpz[i], rhos->x.mpz[L->n-1]));
    }
    
    SPEX_CHECK(SPEX_mpq_init(det2));
    // b2 = L' \ b2
    OK(SPEX_Chol_ltsolve(L, b2));
    
    // x = b2/det 
    mpq_set_num(det2, rhos->x.mpz[L->n-1]);
    
    SPEX_CHECK (SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPQ, b2->m, b->n,
        0, false, true, option)) ;
        
    nz = SPEX_matrix_nnz(b, NULL);
    
    
    for (i = 0; i < nz; i++)
    {
        OK ( SPEX_mpq_set_num( x->x.mpq[i], b2->x.mpz[i]));
        OK ( SPEX_mpq_div( x->x.mpq[i], x->x.mpq[i], det2));
    }
    
    // Permute x
    SPEX_CHECK (SPEX_matrix_allocate(&x2, SPEX_DENSE, SPEX_MPQ, x->m, x->n,
        0, false, true, option)) ;
    
    for (int32_t i = 0; i < x->m; i++)
    {
        for (int32_t j = 0; j < x->n; j++)
        {
            SPEX_CHECK(SPEX_mpq_set( SPEX_2D(x2, S->q[i], j, mpq), SPEX_2D(x, i, j, mpq)));
        }
    }

    // Check solution
    bool check = option->check;
    if (check)
    {
        SPEX_CHECK (SPEX_check_solution (A_orig, x2, b, option)) ;
    }
    
    
    //--------------------------------------------------------------------------
    // Scale the solution if necessary.
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_mpq_init(scale));

    // set the scaling factor scale = A->scale / b->scale
    SPEX_CHECK( SPEX_mpq_div(scale, A->scale, b->scale));

    // Determine if the scaling factor is 1
    int r;
    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, scale, 1, 1));
    nz = x->m * x->n;
    if (r != 0 )
    {
        for (i = 0; i < nz; i++)
        {
            SPEX_CHECK(SPEX_mpq_mul(x2->x.mpq[i], x2->x.mpq[i], scale));
        }
    }
    
    
    // Free memory
    (*x_handle) = x2;
    FREE_WORKSPACE;
    return SPEX_OK;
}
