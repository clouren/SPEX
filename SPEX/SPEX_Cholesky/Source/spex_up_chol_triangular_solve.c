//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Up_Chol_triangular_solve: Sparse symmetric REF Triangular solve
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

# define SPEX_FREE_ALL            \
{                                 \
    SPEX_matrix_free(&x, NULL);   \
    SPEX_FREE(xi);                \
}

#include "spex_chol_internal.h"

/* Purpose: This function performs the symmetric sparse REF triangular solve for
 * the up-looking Cholesky factorization. i,e, LD x = A(1:k-1, k). At the end of
 * this function, the vector x contains the values of the kth row of the integer-
 * preserving matrix L. 
 * 
 * Input arguments of the function:
 *
 * top_output:      A pointer to the beginning of the nonzero pattern. Undefined
 *                  on input, on output xi[top_output..n] contains the beginning
 *                  of the nonzero pattern.
 * 
 * L:               Lower triangular matrix.
 * 
 * A:               Input matrix
 * 
 * k:               Current iteration of the algorithm
 * 
 * xi:              Nonzero pattern. Undefined on input, on output contains teh 
 *                  nonzero pattern of the kth row of L
 * 
 * parent:          Elimination tree
 * 
 * c:               Column pointers of L but they don't point to the top
 *                  position of each column of L. Instead they point to the
 *                  position on each column where the next value of L will be
 *                  grabbed, since at iteration k we need to grab the kth of L
 *                  in order to not recompute those values.
 * 
 * rhos:            Pivot matrix
 * 
 * h:               History vector
 * 
 * x:               Solution of linear system. Undefined on input, on output
 *                  contains the kth row of L.
 */

// Comparison function used for the quicksort in the factorization
// Each iteration of the triangular solve requires that the nonzero pattern
// is sorted prior to numeric operations. This is the helper function for
// c's default qsort
static inline int compare (const void * a, const void * b)
{
    return ( *(int64_t*)a - *(int64_t*)b ) ;
}


SPEX_info spex_up_chol_triangular_solve
(
    //Output
    int64_t* top_output,               // On input NULL. On output contains the
                                       // beginning of nonzero pattern
                                       // The nonzero pattern is contained in 
                                       // xi[top_output...n-1]
    int64_t* xi,                       // Nonzero pattern vector
    SPEX_matrix* x,                    // Solution of system ==> kth row of L
    // Input
    const SPEX_matrix* L,              // Partial L matrix
    const SPEX_matrix* A,              // Input matrix
    const int64_t k,                   // Iteration of algorithm
    const int64_t* parent,             // Elimination tree
    int64_t* c,                        // Column pointers
    const SPEX_matrix* rhos,           // sequence of pivots
    int64_t* h                         // History vector
)
{
    SPEX_info info;
    

    // check inputs
    ASSERT(L->type == SPEX_MPZ);
    ASSERT(L->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(x->type == SPEX_MPZ);
    ASSERT(x->kind == SPEX_DENSE);
    ASSERT(rhos->type == SPEX_MPZ);
    ASSERT(rhos->kind == SPEX_CSC);
    
    int64_t j, i, inew, p, m, top, n = A->n, col;
    int sgn;
    
    ASSERT(n >= 0);
    
    //--------------------------------------------------------------------------
    // Initialize REF Triangular Solve by getting the nonzero patern of x &&
    // obtaining A(:,k)
    //--------------------------------------------------------------------------
    // Obtain the nonzero pattern of the kth row of L by analyzing the
    // elimination tree of A. The indices of these nonzeros are stored in
    // xi[top..n-1]
    SPEX_CHECK(spex_chol_ereach(&top, xi, A, k, parent, c));
    
    // Sort the nonzero pattern using quicksort (reuqired by IPGE unlike in GE)
    qsort(&xi[top], n-top, sizeof(int64_t*), compare); 
        
    // Reset x[i] = 0 for all i in nonzero pattern xi [top..n-1]
    for (i = top; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpz_set_ui(x->x.mpz[xi[i]],0));
    }
    
    // TODO: Add a nonzero diagonal check to the symmertic check (due to the sparse format, this has to be done with a workspace 'marking' array at the same time of the numeric symmerty check)
    // Reset value of x[k]. If the matrix is nonsingular, x[k] will
    // be a part of the nonzero pattern and reset in the above loop.
    // However, in some rare cases, the matrix can be singular but x[k]
    // will be nonzero from a previous iteration. Thus, here we reset
    // x[k] to account for this extremely rare case.
    // Note that this can only happen if the `expert' user did not performed the
    // symmerty/non-zero-diagonal check.
    SPEX_CHECK(SPEX_mpz_set_ui(x->x.mpz[k],0));

    // Reset h[i] = -1 for all i in nonzero pattern
    for (i = top; i < n; i++)
    {
        h[xi[i]] = -1;
    }
    
    // Set x = A(:,k)
    // Note: The if is needed since the columns of A are allowed to be unsorted.
    for (i = A->p[k]; i < A->p[k+1]; i++)
    {
        if (A->i[i] <= k)
        {
            SPEX_CHECK(SPEX_mpz_set(x->x.mpz[A->i[i]], A->x.mpz[i]));
        }
    }
    
    //--------------------------------------------------------------------------
    // Perform the REF Triangular Solve. Note that, unlike the left-looking 
    // Cholesky triangular solve where L is lower trapezoidal, the up-looking 
    // L matrix is actually lower triangular; thus this is a true triangular solve.
    //--------------------------------------------------------------------------
    for (p = top; p < n; p++)
    {   
        // Obtain the index of the current nonzero
        j = xi[p];
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x.mpz[j]));
        if (sgn == 0) continue;    // If x[j] == 0 no work must be done
                
        // Initial history update to finalize x[j] if necessary
        if (h[j] < j-1)
        {
            // History update x[j]: x[j] = x[j]*rhos[j-1]/rhos[h[j]]
            // x[j] = x[j]*rhos[j-1]
            SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[j], x->x.mpz[j], rhos->x.mpz[j-1]));
            if (h[j] > -1)
            {
               // x[j] = x[j] / rhos [ h[j] ]
               SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[j], x->x.mpz[j], 
                                            rhos->x.mpz[h[j]]));
            }
        }
        
        //------------------------------------------------------------------
        // IPGE updates
        //------------------------------------------------------------------
        // ----------- Iterate accross nonzeros in Lij ---------------------
        for (m = L->p[j]+1; m < c[j]; m++)
        {
            i = L->i[m];            // i value of Lij
            if (i > j && i < k)     // Update all dependent x[i] excluding x[k]
            {
                    /*************** If lij==0 then no update******************/
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->x.mpz[m]));
                if (sgn == 0) continue;

                //----------------------------------------------------------
                /************* lij is nonzero, x[i] is zero****************/
                // x[i] = 0 then only perform IPGE update subtraction/division
                //----------------------------------------------------------
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x.mpz[i]));
                if (sgn == 0)
                {
                    // First, get the correct value of x[i] = 0 - lij * x[j]
                    // TODO: Please compare this single fgunction call to doing via two calls: x[i] = lij * x[j], and then flipping the sign of x[i] (flip sign with flip sign function, not multiply by -1) ----all of this check for the case where x[i] is zero
                    // ^ use a list of humongous numbers
                    SPEX_CHECK(SPEX_mpz_submul(x->x.mpz[i], L->x.mpz[m],
                                                 x->x.mpz[j]));
                    // Do a division by the pivot if necessary.
                    if (j >= 1)
                    {
                        // x[i] = x[i] / rho[j-1]
                        SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[i], x->x.mpz[i], 
                                                        rhos->x.mpz[j-1]));
                    }
                    // Update the history value of x[i]
                    h[i] = j;
                    
                }

                //----------------------------------------------------------
                /************ Both lij and x[i] are nonzero****************/
                // x[i] != 0 --> History & IPGE update on x[i]
                //----------------------------------------------------------
                else
                {
                    // TODO for TIM: What do you prefer the block below (longer) or the one commented below (more compact)?
                    // There is no previous pivot
                    if (j < 1)
                    {
                        // History update x[i] = x[i]*rhos[0]
                        SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[i],x->x.mpz[i],
                                                rhos->x.mpz[0]));
                        // x[i] = x[i] - lij x[j]
                        SPEX_CHECK(SPEX_mpz_submul(x->x.mpz[i], L->x.mpz[m], 
                                                    x->x.mpz[j]));
                        // Update the history value of x[i]
                        h[i] = j;
                    }
                    // There is a previous pivot
                    else
                    {
                        // History update if necessary
                        if (h[i] < j - 1)
                        {
                            // x[i] = x[i] * rhos[j-1]
                            SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[i],x->x.mpz[i],
                                                     rhos->x.mpz[j-1]));
                            // Divide by the history pivot only if the history 
                            // pivot is not the rho[-1] (which equals 1) (rho[0] 
                            // in the 1-based logic of othe IPGE algorithm)
                            if (h[i] > -1)
                            {
                                // x[i] = x[i] / rho[h[i]]
                                SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[i],
                                            x->x.mpz[i],rhos->x.mpz[h[i]]));
                            }
                        }
                        // ---- IPGE Update x[i] = (x[i]*rhos[j] - lij*xj) / rho[j-1] ------
                        // x[i] = x[i]*rhos[j]
                        SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[i],x->x.mpz[i],
                                                rhos->x.mpz[j]));
                        // x[i] = x[i] - lij*xj
                        SPEX_CHECK(SPEX_mpz_submul(x->x.mpz[i], L->x.mpz[m], 
                                                    x->x.mpz[j]));
                        // x[i] = x[i] / rho[j-1] 
                        SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[i],x->x.mpz[i],
                                                        rhos->x.mpz[j-1]));
                        // Entry is up to date; 
                        h[i] = j;
                    }
                                        
/*
                    // History update if necessary; That is, if there is the previous pivot is not rho[-1] (which equals 1) (j>0), and the entry is not up-to-date (h[i] < j - 1)
                    if (j > 0 && h[i] < j - 1)
                    {
                        // x[i] = x[i] * rhos[j-1]
                        SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[i], x->x.mpz[i], rhos->x.mpz[j-1]));
                        // Divide by the history pivot only if the history pivot is not the rho[-1] (which equals 1) (rho[0] in the 1-based logic of othe IPGE algorithm)
                        if (h[i] > -1)
                        {
                            // x[i] = x[i] / rho[h[i]]
                            SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[i], x->x.mpz[i], rhos->x.mpz[h[i]]));
                        }
                    }
                    // ---- IPGE Update x[i] = (x[i]*rhos[j] - lij*xj) / rho[j-1] ------
                    // x[i] = x[i]*rhos[0]
                    SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[i], x->x.mpz[i], rhos->x.mpz[j]));
                    // x[i] = x[i] - lij*xj
                    SPEX_CHECK(SPEX_mpz_submul(x->x.mpz[i], L->x.mpz[m], x->x.mpz[j]));
                    // x[i] = x[i] / rho[j-1] 
                    // Only divide by previous pivot if the previous pivot is not 1 (which is always the case in the first IPGE iteration)
                    if (j > 0)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[i], x->x.mpz[i], rhos->x.mpz[j-1]));
                    }
                    // Entry is up to date; 
                    h[i] = j;
 */                   
                }
            }
        }
        // ------ History Update x[k] if necessary -----
        if (h[k] < j - 1)
        {
            // x[k] = x[k] * rho[j-1]
            SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[k],x->x.mpz[k],rhos->x.mpz[j-1]));
            // Divide by the history pivot only if the history pivot is not the 
            // rho[-1] (which equals 1) (rho[0] in the 1-based logic of the 
            // IPGE algorithm)
            if (h[k] > -1)
            {
                // x[k] = x[k] / rho[h[k]]
                SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[k],x->x.mpz[k],
                                              rhos->x.mpz[h[k]]));
            }
        }
        // ---- IPGE Update x[k] = (x[k]*rhos[j] - xj*xj) / rho[j-1] ------
        // x[k] = x[k] * rho[j]
        SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[k],x->x.mpz[k],rhos->x.mpz[j]));
        // x[k] = x[k] - xj*xj
        SPEX_CHECK(SPEX_mpz_submul(x->x.mpz[k], x->x.mpz[j], x->x.mpz[j]));
        // Only divide by previous pivot if the previous pivot is not 1 (which 
        // is always the case in the first IPGE iteration)
        if (j > 0)
            // x[k] = x[k] / rho[j-1] 
            SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[k],x->x.mpz[k],
                                            rhos->x.mpz[j-1])); 
        // Entry is up to date; 
        h[k] = j;   
    }
    //----------------------------------------------------------
    // At this point, x[k] has been updated throughout the 
    // triangular solve. The last step is to make sure x[k]
    // has its correct final value. Thus, a final history
    // update is done to x[k] if necessary
    //----------------------------------------------------------
    if (h[k] < k-1)
    {
        // x[k] = x[k] * rhos[k-1]
        SPEX_CHECK(SPEX_mpz_mul(x->x.mpz[k], x->x.mpz[k], rhos->x.mpz[k-1]));
        // Only divide by previous pivot if the previous pivot is not 1 (which 
        // is always the case in the first IPGE iteration)
        if (h[k] > -1)
        {
            // x[k] = x[k] / rhos[h[k]]
            SPEX_CHECK(SPEX_mpz_divexact(x->x.mpz[k], x->x.mpz[k], 
                                           rhos->x.mpz[ h[k]]));
        }
    }
    // Output the top of the nonzero pattern
    *top_output = top;
    return SPEX_OK;
}