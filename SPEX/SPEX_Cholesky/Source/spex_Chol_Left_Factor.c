//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_Left_Factor: Left-looking REF Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define FREE_WORKSPACE              \
    SPEX_matrix_free(&x, NULL);     \
    SPEX_FREE(c);                   \
    SPEX_FREE(xi);                  \
    SPEX_FREE(h);                   \
    SPEX_FREE(post);                \


#include "spex_chol_internal.h"  

/* Purpose: This function performs the left-looking REF Cholesky factorization.
  * In order to compute the L matrix, it performs n iterations of a sparse REF symmetric
 * triangular solve function which, at each iteration, computes the kth column of L. 
 * 
 * Note that, importantly, this function assumes that A has been permuted
 * prior to factorization.
 * 
 * Input arguments:
 * 
 * L_handle:    A handle to the L matrix. Null on input. On output, contains a pointer to the 
 *              L matrix
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. On output,
 *              contains the elimination tree and estimate number of nonzeros in L.
 * 
 * rhos_handle: A handle to the sequence of pivots. NULL on input. On output, contains a pointer
 *              to the pivots matrix.
 *
 * A:           The user's permuted input matrix
 * 
 * option:      Command options
 * 
 */
SPEX_info spex_Chol_Left_Factor      
(
    // Output
    SPEX_matrix** L_handle,     // Lower triangular matrix. NULL on input
    SPEX_matrix** rhos_handle, // Sequence of pivots. NULL on input.
    SPEX_Chol_analysis* S,     // Symbolic analysis struct that contains elimination tree of A, column pointers of L, 
                                //exact number of nonzeros of L and permutation used
    // Input
    const SPEX_matrix* A,      // Matrix to be factored   
    const SPEX_options* option //command options
)
{
    SPEX_info info;
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!L_handle || !S || !rhos_handle || !option )
    {
        return SPEX_INCORRECT_INPUT;
    }
    ASSERT(*L_handle==NULL);
    ASSERT(*rhos_handle==NULL);

    int64_t anz;
    SPEX_CHECK(SPEX_matrix_nnz (&anz, A, option)) ;
    
    if (anz < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*rhos_handle) = NULL ;
        
    //--------------------------------------------------------------------------
    // Declare and initialize workspace 
    //--------------------------------------------------------------------------
    
    SPEX_matrix *L = NULL ;
    SPEX_matrix *rhos = NULL ;
    int64_t *xi = NULL ;
    int64_t *h = NULL ;
    SPEX_matrix *x = NULL ;

    // Declare variables
    int64_t n = A->n, top, i, j, col, loc, lnz = 0, unz = 0, jnew, k;
    size_t size;

    // Post and c are used to compute the elimination tree
    int64_t* post = NULL;
    int64_t* c = NULL;


    
    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed into the triangular
    // solve algorithm
    h = (int64_t*) SPEX_malloc(n* sizeof(int64_t));

    // xi is the global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L and U
    // for the triangular solve.
    xi = (int64_t*) SPEX_malloc(2*n* sizeof(int64_t));

    if (!h || !xi)
    {
        FREE_WORKSPACE;
        return SPEX_OUT_OF_MEMORY;
    }
    // initialize workspace history array
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
    }

        
    // Obtain the elimination tree of A
    SPEX_CHECK(spex_Chol_etree(&S->parent, A));
    
    // Postorder the elimination tree of A
    SPEX_CHECK( spex_Chol_post(&post, S->parent, n));    // Postorder the tree
    
    // Get the column counts of A
    SPEX_CHECK( spex_Chol_counts(&c, A, S->parent, post));
    
    // Set the column pointers of L
    S->cp = (int64_t*) SPEX_malloc( (n+1)*sizeof(int64_t*));
    SPEX_CHECK( SPEX_cumsum(S->cp, c, n));
   
    S->lnz = S->cp[n];
   
    //--------------------------------------------------------------------------
    // allocate and initialize the workspace x
    //--------------------------------------------------------------------------

    // SPEX utilizes arbitrary sized integers which can grow beyond the
    // default 64 bits allocated by GMP. If the integers frequently grow, GMP
    // can get bogged down by performing intermediate reallocations. Instead,
    // we utilize a larger estimate on the workspace x vector so that computing
    // the values in L and U do not require too many extra intermediate calls to
    // realloc.
    //
    // Note that the estimate presented here is not an upper bound nor a lower
    // bound.  It is still possible that more bits will be required which is
    // correctly handled internally.
    int64_t estimate = 64 * SPEX_MAX (2, ceil (log2 ((double) n))) ;
    
    // Create x, a global dense mpz_t matrix of dimension n*1. Unlike rhos, the
    // second boolean parameter is set to false to avoid initializing
    // each mpz entry of x with default size.  It is initialized below.
    SPEX_CHECK (SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, /* do not initialize the entries of x: */ false, option));
    
    // Create rhos, a global dense mpz_t matrix of dimension n*1. 
    SPEX_CHECK (SPEX_matrix_allocate(&rhos, SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));
    
    if (!x || !rhos)
    {
        FREE_WORKSPACE;
        return SPEX_OUT_OF_MEMORY;
    }
    
    // initialize the entries of x
    for (i = 0; i < n; i++)
    {
        // Allocate memory for entries of x
        SPEX_CHECK(SPEX_mpz_init2(x->x.mpz[i], estimate));
    }
    
    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
    
    // Since we are performing a left-looking factorization, we pre-allocate L
    // by performing a symbolic version of the factorization and obtaining the 
    // exact nonzero pattern of L.
    // That said, the inidividual (x) values of L are not allocated. Instead,
    // a more efficient method to allocate these values is done inside the 
    // factorization to reduce memory usage.
    
    SPEX_CHECK(spex_Chol_Pre_Left_Factor(&L, xi, A, S->parent, S, c));
        
    // Set the column pointers of L
    for (k = 0; k < n; k++) 
        L->p[k] = c[k] = S->cp[k];
    

    //--------------------------------------------------------------------------
    // Perform the factorization
    //--------------------------------------------------------------------------
        
    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (1:n in standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        SPEX_CHECK (spex_Left_Chol_triangular_solve(&top, x, xi, L, A, k, rhos, h, S->parent, c));

        // Set the pivot element If this element is equal to zero, no pivot element exists and
        // the matrix is either not SPD or singular
        if (mpz_sgn(x->x.mpz[k]) != 0)
        {
            SPEX_CHECK(SPEX_mpz_set(rhos->x.mpz[k], x->x.mpz[k]));
        }
        else
        {
            FREE_WORKSPACE;
            return SPEX_SINGULAR;
        }
        //----------------------------------------------------------------------
        // Add the nonzeros to the L matrix
        //----------------------------------------------------------------------
        for (j = top; j < n; j++)
        {
            // Index of x[i]
            jnew = xi[j];
            if (jnew >= k)
            {
                // Find the size of x[j]
                size = mpz_sizeinbase(x->x.mpz[jnew],2);
                    
                // GMP manual: Allocated size should be size+2
                SPEX_CHECK(SPEX_mpz_init2(L->x.mpz[lnz], size+2));
                
                // Place the x value of this nonzero in row jnew
                SPEX_CHECK(SPEX_mpz_set(L->x.mpz[lnz],x->x.mpz[jnew]));
                // Increment lnz
                lnz += 1;
            }
        }
    }
    // Finalize L->p
    L->p[n] = S->lnz;
       
    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    (*L_handle) = L;
    (*rhos_handle) = rhos;
    FREE_WORKSPACE;
    return SPEX_OK;
}
