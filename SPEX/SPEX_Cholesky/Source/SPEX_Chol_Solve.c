//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_Solve: Solve the SPD linear system after factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#define SPEX_FREE_WORKSPACE        \
    SPEX_matrix_free(&b2, option); \
    SPEX_matrix_free(&x, option); \
    SPEX_MPQ_CLEAR(scale);         

# define SPEX_FREE_ALLOCATION      \
    SPEX_FREE_WORKSPACE            \
    SPEX_matrix_free(&x, NULL);    

#include "spex_chol_internal.h"
    
/* Purpose: This function solves the linear system LDL' x = b.
 *
 * Input arguments:
 * 
 * x_handle:        A handle to the solution matrix. On input this is NULL,
 *                  on output x_handle contains a pointer to the solution
 *                  vector(s)
 * 
 * PAP:             Permuted version of the input matrix
 * 
 * A:               Nonpermuted (original) input matrix
 * 
 * b:               Right hand side vector(s)
 * 
 * rhos:            Sequence of factorization pivots
 * 
 * L:               Lower triangular matrix
 * 
 * S:               Symbolic analysis struct (contains row/column permutation
 *                  and its inverse
 * 
 * option:          Command options
 * 
 */

SPEX_info SPEX_Chol_Solve       // solves the linear system LDL' x = b
(
    // Output
    SPEX_matrix** x_handle,           // On input: undefined.
                                      // On output: Rational solution (SPEX_MPQ)
                                      // to the system. 
    // Input
    //const SPEX_matrix* PAP,           // Input matrix (permuted)
    //const SPEX_matrix* A,             // Input matrix (unpermuted)
    const SPEX_factorization* F,      // Cholesky factorization
    const SPEX_matrix* b,             // Right hand side vector
    //const SPEX_Chol_analysis* S,      // Symbolic analysis struct. contains
                                      // elimination tree of A,  
                                      // column pointers of L, exact number of
                                      // nonzeros of L and permutation used
    const SPEX_options* option        // command options
)
{
    SPEX_info info;

    if (!spex_initialized()) return SPEX_PANIC;
  
    // Check the inputs
    ASSERT(!x_handle);
    ASSERT(b->type == SPEX_MPZ);
    ASSERT(b->kind == SPEX_DENSE);

    //TOASK do we want to check the type and kind of things inside F??
    if (!x_handle || b->type != SPEX_MPZ || b->kind != SPEX_DENSE)
    {
        return SPEX_INCORRECT_INPUT;
    }
    
    int64_t i, j, k, n = F->L->n, nz;
    
    // Scale is the scaling factor for the solution vectors.
    // When the forward/backsolve is complete, the entries in
    // x are rational, but are solving the scaled linear system
    // A' x = b' (that is if A had input which was rational or floating point
    // and had to be converted to integers). Thus, the scale here is used
    // to convert x into into the actual solution of A x = b. 
    // Mathematically, scale = PAP->scale / b->scale
    mpq_t scale;
    SPEX_MPQ_SET_NULL(scale);
  
    // det is the determinant of the PAP matrix. It is obtained for free
    // from the SPEX Cholesky factorization det = rhos[n-1] = L[n,n]
    mpz_t* det = NULL;

    // Declare workspace and output
    // x is the unpermuted solution vector
    SPEX_matrix *x = NULL;
    // x2 is the permuted final solution vector returned to the user
    SPEX_matrix *x2 = NULL;
    // b2 is the permuted right hand side vector(s)
    SPEX_matrix *b2 = NULL;
    //TOASK left_lu_permute_b ???
    // Allocate memory for b2
    SPEX_CHECK(SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPZ, b->m, b->n, 
                                        b->m*b->n, false, true, option));
    // Set b2[i,j] = b[pinv[i],j] to account for permutations applied to A
    for (i = 0; i < b->m; i++)
    {
        for (j = 0; j < b->n; j++)   
        {
            SPEX_CHECK(SPEX_mpz_set(SPEX_2D(b2, F->Pinv_perm[i], j, mpz),
                              SPEX_2D(b, i, j, mpz)));
        }
    }
    // Forward substitution, b2 = L \ b2. Note that b2 is overwritten    
    SPEX_CHECK(spex_Chol_forward_sub(b2, F->L, F->rhos));
    // Set the value of the determinant det = rhos[n-1] 
    det = &(F->rhos->x.mpz[F->L->n-1]);
    size_t bitsDet = mpz_sizeinbase((*det),2); //CLUSTER
    //gmp_printf("%llu, ", (bitsDet)); //CLUSTER
    
    //TOASK SPEX_matrix_mul vs SPEX_mpz_mul + for                          
    // Apply the determinant to b2, b2 = det*b2
    SPEX_CHECK(SPEX_matrix_nnz(&nz, b, NULL));
    for (i = 0; i < nz; i++)
    {
        SPEX_CHECK(SPEX_mpz_mul(b2->x.mpz[i], b2->x.mpz[i], (*det) ));
    }
    
    // Backsolve, b2 = L' \ b2. Note that, again, b2 is overwritten
    SPEX_CHECK(spex_Chol_backward_sub(b2, F->L));

    // Allocate memory for the (real) solution vector x = b2/det
    SPEX_CHECK(SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPQ, b2->m, b->n,
        0, false, true, option)) ;
/*
    size_t bitsSmall=mpz_sizeinbase(b2->x.mpz[0],2); //CLUSTER
    size_t bitsLarge=0; //CLUSTER
    size_t bitsTemp; //CLUSTER
    int64_t bitsSum=0; //CLUSTER
    // Set x[i] = b2[i]/det
*/
    for (i = 0; i < nz; i++)
    {
        /*
        bitsTemp=mpz_sizeinbase(b2->x.mpz[i],2); //CLUSTER
        if(bitsTemp>bitsLarge) //CLUSTER
        { 
            bitsLarge=bitsTemp; //CLUSTER
        }
        else if(bitsTemp<bitsSmall) //CLUSTER
        {
                bitsSmall=bitsTemp; //CLUSTER
        }
        bitsSum=bitsSum+(int64_t)bitsTemp; //CLUSTER
*/
        SPEX_CHECK(SPEX_mpq_set_num(x->x.mpq[i], b2->x.mpz[i]));
        /*mpq_t temp;
        mpq_init(temp);
        mpq_set_z(temp,(*det));
        mpq_div(x->x.mpq[i],x->x.mpq[i],temp);*/ //THIS WORKS

        //SPEX_CHECK(SPEX_mpq_set_den(x->x.mpq[i], (*det) ));
        SPEX_CHECK(SPEX_mpz_set(mpq_denref(x->x.mpq[i]), (*det)));
        SPEX_CHECK(SPEX_mpq_canonicalize(x->x.mpq[i])); //remove all common factors 
    }
    /*
    gmp_printf("%llu, %llu, %llu, ", bitsSmall, bitsLarge, bitsSum); //CLUSTER
    */
    // Allocate memory for x2 which is the permuted version of x
    SPEX_CHECK(SPEX_matrix_allocate(&x2, SPEX_DENSE, SPEX_MPQ, x->m, x->n,
                                    0, false, true, option));
    
    // Apply the row permutation x2[ p[i], j] = x[i,j]
    for (i = 0; i < x->m; i++)
    {
        for (j = 0; j < x->n; j++)
        {
            SPEX_CHECK(SPEX_mpq_set(SPEX_2D(x2, F->Q_perm[i], j, mpq), SPEX_2D(x, i,
                                                                    j, mpq)));
        }
    }


    // Check solution
    // TODO: Shouldnt this be removed if/when Jinhao changes function?
    //TODO where are we adding this??
    /*if (option->check)
    {
        SPEX_CHECK(SPEX_check_solution(A, x2, b, option));
    }
    */
    
    //--------------------------------------------------------------------------
    // Scale the solution if necessary.
    //--------------------------------------------------------------------------

    //TODO doublecheck scaling process
    SPEX_CHECK(SPEX_mpq_init(scale));

    // set the scaling factor scale = PAP->scale / b->scale
    //TODO TOASK double check PAP scale is the same thing as F->scale for A
    SPEX_CHECK(SPEX_mpq_div(scale, F->scale_for_A, b->scale));

    // Apply scaling factor, but ONLY if it is different from a scaling factor
    // to 1
    int r;
    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, scale, 1, 1));
    if (r != 0)
    {
        nz = x->m * x->n;
        for (i = 0; i < nz; i++)
        {
            SPEX_CHECK(SPEX_mpq_mul(x2->x.mpq[i], x2->x.mpq[i], scale));
        }
    }
    
    // Set output, free memory
    (*x_handle) = x2;
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}
