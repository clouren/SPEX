//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_Up_Factor: Up-looking REF Cholesky factorization
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

/* Purpose: This function performs the up-looking REF Cholesky factorization.
 * In order to compute the L matrix, it performs n iterations of a sparse
 * REF triangular solve function. The up-looking algorithm computes L one row
 * at a time.
 * 
 * Note that, importantly, it is expected that the input matrix A has already
 * been permuted. That is, the L matrix computed here is directly the factorization
 * of the input A.
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
 * left:        A boolean parameter which tells the function whether it is performing a left-looking
 *              or up-looking factorization. If this bool is true, a left-looking factorization
 *              is done, otherwise the up-looking factorization is done.
 * 
 * option:      Command options
 * 
 */

SPEX_info spex_Chol_Up_Factor      
(
    // Output
    SPEX_matrix** L_handle,    // Lower triangular matrix. NULL on input
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
    
    // A must be CSC and MPZ
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

    // Post and c are vectors utilized for the construction of the elimination
    // tree
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
    
    // Postorder the tree of A
    SPEX_CHECK( spex_Chol_post(&post, S->parent, n));
    
    // Get the column counts of A
    SPEX_CHECK( spex_Chol_counts(&c, A, S->parent, post));
    
    // Get the column pointers of L
    S->cp = (int64_t*) SPEX_malloc( (n+1)*sizeof(int64_t*));
    SPEX_CHECK( SPEX_cumsum(S->cp, c, n));
   
    // Set the number of nonzeros in L
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
    // The bound given in the paper is that the number of bits is <= n log sigma
    // where sigma is the largest entry in A. Instead of using this bound, we use
    // a very rough estimate of 64*max(2, log (n))
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
    
    // Since we are performing an up-looking factorization, we allocate
    // L without initializing each entry.
    // Note that, the inidividual (x) values of L are not allocated. Instead,
    // a more efficient method to allocate these values is done inside the 
    // factorization to reduce memory usage.
    
    SPEX_CHECK (SPEX_matrix_allocate(&L, SPEX_CSC, SPEX_MPZ, n, n, S->lnz,
            false, false, option));
    
    // Set the column pointers of L
    for (k = 0; k < n; k++) 
        L->p[k] = c[k] = S->cp[k];
    

    //--------------------------------------------------------------------------
    // Perform the up-looking factorization
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (1:n in standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        SPEX_CHECK ( spex_Up_Chol_triangular_solve(&top, xi, x, L, A, k, S->parent, c, rhos, h));
        
        // If x[k] is nonzero that is the pivot. if x[k] == 0 then matrix is 
        // either not SPD or singular.
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
        // Add the nonzeros to L
        //----------------------------------------------------------------------
        int64_t p = 0;
        for (j = top; j < n; j++)
        {
            // Obtain the column index of x[j]
            jnew = xi[j];
            if (jnew == k) continue;
            
            // Determine the column that x[j] belongs in
            p = c[jnew]++;
            
            // Place the i location of this nonzero. This should always be k
            // because at iteration k, the up-looking algorithm computes row k
            // of L
            L->i[p] = k;
            
            // Find the number of bits of x[j]
            size = mpz_sizeinbase(x->x.mpz[jnew],2);
            
            // GMP manual: Allocated size should be size+2
            SPEX_CHECK(SPEX_mpz_init2(L->x.mpz[p], size+2));
            
            // Place the x value of this nonzero
            SPEX_CHECK(SPEX_mpz_set(L->x.mpz[p],x->x.mpz[jnew]));
        }
        // Now, place L(k,k)
        p = c[k]++;
        L->i[p] = k;
        size = mpz_sizeinbase(x->x.mpz[k], 2);
        SPEX_CHECK(SPEX_mpz_init2(L->x.mpz[p], size+2));
        SPEX_CHECK(SPEX_mpz_set(L->x.mpz[p], x->x.mpz[k]));
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
