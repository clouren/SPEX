//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_permute_A: Symmetric permutation of matrix A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/*#define SPEX_FREE_ALL             \
{                                 \
    SPEX_matrix_free(pinv, NULL); \
}*/

#include "spex_chol_internal.h"

/* Purpose: Given the row/column permutation P stored in S, permute the matrix
 * A and return PAP'
 * Input arguments:
 * 
 * PAP_handle:   The user's permuted input matrix. 
 * 
 * A:            The user's input matrix
 * 
 * S:            Symbolic analysis struct for Cholesky factorization. 
 *               Contains row/column permutation of A

 FIXME: why pass in S?  Why not pass in both the permutation P and
 its inverse?  If the inverse passed in is a NULL pointer, then
 compute it, use it, then discard it.  Also, why are both permutations
 Q_perm and Pinv_perm needed?  Only Pinv_perm is required.  See cs_symperm in
 CSparse.  Finally, we should use P_perm not Q_perm.

 * 
 */
SPEX_info spex_chol_permute_A
(
    //Output
    SPEX_matrix** PAP_handle,  // On input: undefined
                               // On output: contains the permuted matrix
    //Input
    const SPEX_matrix* A,      // Input matrix
    const bool numeric,        // True if user wants to permute pattern and 
                               // numbers, false if only pattern
    
    //Input/Ouput
    SPEX_symbolic_analysis* S      // Symbolic analysis struct that contains 
                               // row/column permutations
)
{
    SPEX_info info;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    ASSERT(A != NULL);
    ASSERT(S != NULL);
    ASSERT(PAP_handle != NULL);
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->kind == SPEX_CSC);

    if (!PAP_handle || !S || !A || A->type != SPEX_MPZ || A->kind != SPEX_CSC)
    {
        return SPEX_INCORRECT_INPUT;
    }     

    // Create indices and pinv, the inverse row permutation
    int64_t j, k, t, index, nz = 0, n = A->n;
    //int64_t* pinv = NULL;

    // FIXME: this should not be here; S should be read-only
    // move to spex_chol_preorder
    // Allocate pinv
    if (!(S->Pinv_perm))
    {
        S->Pinv_perm = (int64_t*)SPEX_calloc(n, sizeof(int64_t));
        if(!(S->Pinv_perm))
        {
            return SPEX_OUT_OF_MEMORY;
        }
    }
    // Populate pinv
    for (k = 0; k < n; k++)
    {
        // FIXME: this should not be here; S should be read-only,
        // and S->Pinv_perm should only be computed once.  If this function
        // is called twice, it recomputes S->Pinv_perm every time.
        index = S->Q_perm[k];
        S->Pinv_perm[index] = k;
    }


    //TODO maybe there is a more elegant way? this happens because otherwise we never free the Pinv created in analysis
    /*if((S->Pinv_perm)!=NULL)
    { 
        SPEX_FREE(S->Pinv_perm);
    }*/ //using this makes example simple fail because of memory error
    //S->Pinv_perm=pinv; //TODO should be a full copy

    // Allocate memory for PAP which is a permuted copy of A
    SPEX_matrix* PAP = NULL;
    //SPEX_CHECK(SPEX_matrix_allocate(&PAP, SPEX_CSC, SPEX_MPZ, n, n, A->p[n], false, true, NULL));
    SPEX_CHECK(SPEX_matrix_allocate(&PAP, SPEX_CSC, SPEX_MPZ, n, n, A->p[n], true, false, NULL));
    PAP->p=(int64_t*)SPEX_malloc((n+1)*sizeof(int64_t));
    PAP->i=(int64_t*)SPEX_malloc((A->p[n])*sizeof(int64_t));
    PAP->p_shallow = false ;
    PAP->i_shallow = false ; //TODO FIXME still feels like patchwork, figure out why

    // FIXME: PAP->x.mpz is a different kind of shallow
    
    if(numeric)
    {
       PAP->x.mpz=(mpz_t*)SPEX_malloc((A->p[n])*sizeof(mpz_t));

        // Set PAP scale
        SPEX_CHECK(SPEX_mpq_set(PAP->scale, A->scale));
    
        // Populate the entries in PAP
        for (k = 0; k < n; k++)
        {
            // Set the number of nonzeros in the kth column of PAP
            PAP->p[k] = nz;
            // FIXME: old comment here, out of date:
            // Column k of PAP is equal to column S->p[k] of A. j is the starting
            // point for nonzeros and indices for column S->p[k] of A
            j = S->Q_perm[k];
            // Iterate across the nonzeros in column S->p[k]
            for (t = A->p[j]; t < A->p[j+1]; t++)
            {
                // Set the nonzero value and location of the entries in column k of PAP
                // NOTE: this is shallow.   Provide option for a deep copy of the values?
                // FIXME: call an mpz_* function to do the pointer assignment?
                (*(PAP->x.mpz[nz]))=(*(A->x.mpz[t])); 
                //SPEX_CHECK(SPEX_mpz_set(PAP->x.mpz[nz], A->x.mpz[t]));
                // Row i of this nonzero is equal to pinv[A->i[t]]
                PAP->i[nz] = S->Pinv_perm[ A->i[t] ];
                // Move to the next nonzero element of PAP
                nz++;
            }
        }
    }
    else
    {

        PAP->x.mpz= NULL ;
        PAP->x_shallow = true ;

        // Populate the entries in PAP
        for (k = 0; k < n; k++)
        {
            // Set the number of nonzeros in the kth column of PAP
            PAP->p[k] = nz;
            // Column k of PAP is equal to column S->p[k] of A. j is the starting
            // point for nonzeros and indices for column S->p[k] of A
            j = S->Q_perm[k];
            // Iterate across the nonzeros in column S->p[k]
            for (t = A->p[j]; t < A->p[j+1]; t++)
            {
                // Row i of this nonzero is equal to pinv[A->i[t]]
                PAP->i[nz] = S->Pinv_perm[ A->i[t] ];
                // Move to the next nonzero element of PAP
                nz++;
            }
        }
    } 

    // Finalize the last column of PAP
    PAP->p[n] = nz;

    // Set output, return success
    (*PAP_handle) = PAP;
    return SPEX_OK;
}
