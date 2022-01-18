//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_Up_Factor: Up-looking REF Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE         \
    SPEX_matrix_free(&x, NULL);     \
    SPEX_FREE(c);                   \
    SPEX_FREE(xi);                  \
    SPEX_FREE(h);                   \
    SPEX_FREE(post);                // TOCHECK all the allocations and the free are driving me crazy

//TOASK idk if this free F goes here
# define SPEX_FREE_ALLOCATION        \
    SPEX_FREE_WORKSPACE              \
    //SPEX_factorization_free(&F, option);    //TOCHECK need to declare F before using the free 

#include "spex_chol_internal.h"  

/* Purpose: This function performs the up-looking REF Cholesky factorization.
 * In order to compute the L matrix, it performs n iterations of a sparse REF 
 * symmetric triangular solve function which, at each iteration, computes the
 * kth row of L. 
 * 
 * Importantly, this function assumes that A has already been permuted.
 * 
 * Input arguments of the function:
 * 
 * L_handle:    A handle to the L matrix. Null on input.
 *              On output, contains a pointer to the L matrix.
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. 
 *              On input it contains information that is not used in this 
 *              function such as the row/column permutation
 *              On output it contains the elimination tree and 
 *              the number of nonzeros in L.
 * 
 * rhos_handle: A handle to the sequence of pivots. NULL on input. 
 *              On output it contains a pointer to the pivots matrix.
 *
 * A:           The user's permuted input matrix
 * 
 * option:      Command options
 * 
 */


SPEX_info spex_Chol_Up_Factor      
(
    // Output
    SPEX_matrix** L_handle,    // Lower triangular matrix. NULL on input.
    SPEX_matrix** rhos_handle, // Sequence of pivots. NULL on input.
    //SPEX_factorization **F_handle,
    // Input/Output
    SPEX_symbolic_analysis* S,     // Symbolic analysis struct containing the
                               // elimination tree of A, the column pointers of
                               // L, and the exact number of nonzeros of L.
    // Input
    const SPEX_matrix* A,      // Matrix to be factored   
    const SPEX_options* option // command options
)
{
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    
    if (!L_handle || !S || !option || !A)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    ASSERT(A->type == SPEX_MPZ) ;
    ASSERT(A->kind == SPEX_CSC) ;

    // Check the number of nonzeros in A
    int64_t anz;
    // SPEX enviroment is checked to be init'ed and A is a SPEX_CSC matrix that
    // is not NULL, so SPEX_matrix_nnz must return SPEX_OK
    SPEX_info info = SPEX_matrix_nnz (&anz, A, option) ;
    ASSERT(info == SPEX_OK); //REMINDER in SPEX_CHECK if info!=SPEX_OK it would Free alloc (which is why we can't use SPEX_CHECK here)
        //TODO remove reminder (but not until there is no chance of making same mistake)
    
    if (anz < 0)
    {
        return SPEX_INCORRECT_INPUT ;
    }

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
    int sgn;
    size_t size;

    // Post & c are arrays used for the construction of the elimination tree
    int64_t* post = NULL;
    int64_t* c = NULL;

    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed into the triangular
    // solve algorithm
    h = (int64_t*) SPEX_malloc(n*sizeof(int64_t));

    // xi serves as a global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L
    // for the triangular solve.
    xi = (int64_t*) SPEX_malloc(2*n*sizeof(int64_t));

    if (!h || !xi)
    {
        SPEX_FREE_WORKSPACE;
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
    SPEX_CHECK(spex_Chol_post(&post, S->parent, n));
    
    // Get the column counts of A
    SPEX_CHECK(spex_Chol_counts(&c, A, S->parent, post));
    
    // Get the column pointers of L
    S->cp = (int64_t*) SPEX_malloc((n+1)*sizeof(int64_t*));
    SPEX_CHECK(SPEX_cumsum(S->cp, c, n));
   
    // Set the exact number of nonzeros in L
    S->lnz = S->cp[n];
   
    //--------------------------------------------------------------------------
    // Allocate and initialize the workspace x
    //--------------------------------------------------------------------------

    // SPEX utilizes arbitrary sized integers which can grow beyond the
    // default 64 bits allocated by GMP. If the integers frequently grow, GMP
    // can get bogged down by performing intermediate reallocations. Instead,
    // we utilize a larger estimate on the workspace x vector so that computing
    // the values in L and U do not require too many extra intermediate calls to
    // realloc.
    //
    // The bound given in the paper is that the number of bits is <= n log sigma
    // where sigma is the largest entry in A. Because this bound is extremely 
    // pessimistic, instead of using this bound, we use a very rough estimate:
    // 64*max(2, log (n))
    //
    // Note that the estimate presented here is not an upper bound nor a lower
    // bound.  It is still possible that more bits will be required which is
    // correctly handled internally.
    int64_t estimate = 64 * SPEX_MAX (2, ceil (log2 ((double) n))) ;
    
    // Create x, a "global" dense mpz_t matrix of dimension n*1 (i.e., it is 
    // used as workspace re-used at each iteration). The second boolean
    // parameter is set to false, indicating that the size of each mpz entry
    // will be initialized afterwards (and should not be initialized with the
    // default size)
    SPEX_CHECK (SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, /* do not initialize the entries of x: */ false, option));
    
    // Create rhos, a "global" dense mpz_t matrix of dimension n*1. 
    // As inidicated with the second boolean parameter true, the mpz entries in
    // rhos are initialized to the default size (unlike x).
    SPEX_CHECK (SPEX_matrix_allocate(&(rhos), SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));
    
    if (!x || !rhos)
    {
        SPEX_FREE_WORKSPACE;
        return SPEX_OUT_OF_MEMORY;
    }
    
    // initialize the entries of x
    for (i = 0; i < n; i++)
    {
        // Allocate memory for entries of x to be estimate bits
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
    
    SPEX_CHECK(SPEX_matrix_allocate(&(L), SPEX_CSC, SPEX_MPZ, n, n, S->lnz,
                                    false, false, option));
    
    // Set the column pointers of L
    for (k = 0; k < n; k++) 
    {
        L->p[k] = c[k] = S->cp[k];
    }
    

    //--------------------------------------------------------------------------
    // Perform the up-looking factorization
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (1:n in standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        SPEX_CHECK(spex_Up_Chol_triangular_solve(&top, xi, x, L, A, k, S->parent, 
                                                 c, rhos, h));
  
        // If x[k] is nonzero chose it as pivot. Otherwise, the matrix is 
        // not SPD (indeed, it may even be singular).
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x.mpz[k])); //TODO valgrind says sgn is not initialized, what?
        if (sgn != 0)
        {
            SPEX_CHECK(SPEX_mpz_set(rhos->x.mpz[k], x->x.mpz[k]));
        }
        else
        {
            SPEX_FREE_WORKSPACE;
            // TODO: We need to create the error SPEX_NOTSPD DONE
            // When this is done we also need to change SPEX_Backslash.c 
            // so that it correctly classifies what happens.
            return SPEX_NOTSPD; 
        }
            
        //----------------------------------------------------------------------
        // Add the nonzeros (i.e. x) to L
        //----------------------------------------------------------------------
        int64_t p = 0;
        for (j = top; j < n; j++)
        {
            // Obtain the row index of x[j]
            jnew = xi[j];
            if (jnew == k) continue;
            
            // Determine the column where x[j] belongs to
            p = c[jnew]++;
            
            // Place the i index of this nonzero. This should always be k because 
            // at iteration k, the up-looking algorithm computes row k of L
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
    // Free memory and set output
    //--------------------------------------------------------------------------
    (*L_handle) = L;
    (*rhos_handle) = rhos;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
