//------------------------------------------------------------------------------
// SPEX_Update/spex_update_internal: include file for internal
// use in SPEX_Update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_Update.h instead.

// uncommemd to help debugging
// #define SPEX_DEBUG

#ifndef SPEX_UPDATE_INTERNAL_H
#define SPEX_UPDATE_INTERNAL_H

#include "spex_util_internal.h"
#include "SPEX_Left_LU.h"
#include "SPEX_Update.h"

#ifdef SPEX_DEBUG
// redefine SPEX_CHECK to print the file name and line
#ifdef SPEX_CHECK
#undef SPEX_CHECK
#endif
#define SPEX_CHECK(method)      \
{                               \
    info = (method) ;           \
    if (info != SPEX_OK)        \
    {                           \
        printf("file %s line %d\n",__FILE__,__LINE__);\
        SPEX_FREE_ALL ;         \
        return (info) ;         \
    }                           \
}

#ifndef GOTCHA
#define GOTCHA \
    printf ("GOTCHA: %s, line %d\n", __FILE__, __LINE__);
#endif

#endif

#ifndef GOTCHA
    #define GOTCHA
#endif

// ============================================================================
//                           Internal Functions
// ============================================================================

// define a new type for scattered vector, which is acutally just SPEX_vector
#define spex_scattered_vector SPEX_vector
#define spex_scattered_vector_alloc SPEX_vector_allocate
#define spex_scattered_vector_free SPEX_vector_free

// build scattered vector for given sparse vector and search for inext/jnext if
// requested.
SPEX_info spex_update_get_scattered_v
(
    // output
    spex_scattered_vector **sv_handle,// output vector in scattered form
    int64_t *next,               // the smallest col/row index of non-pivot nz.
                                 // If next == NULL, searching is not performed.
    // input
    SPEX_vector *v,              // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const int64_t k,             // column index of v in L, or row index in U.
                                 // Ignored if next == NULL.
    const int64_t *perm_inv,     // inverse of permutation applied on v.
                                 // This can be NULL if next == NULL.
    const bool keep_v,           // indicate if the mpz values should be kept
    const SPEX_options *option   // command options
);


// perform column permutation pivot update
SPEX_info spex_update_cppu
(
    SPEX_matrix *L,              // matrix L
    SPEX_matrix *U,              // matrix U
    SPEX_matrix *rhos,           // array of scaled pivots
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,              // the index of first off-diag entry in
                                 // column k of L
    int64_t *jnext,              // the index of first off-diag entry in
                                 // row k of U
    int64_t *h,                  // allocated vector that can be used for
                                 // history vector. All entries are maintained
                                 // to be >= -1
    int64_t *Q,                  // column permutation
    int64_t *Q_inv,              // inverse of column permutation
    int64_t *P,                  // row permutation (unchanged on output)
    int64_t *P_inv,              // inverse of row permutation (unchanged
                                 // on output)

    // the col-wise nnz pattern of U can be NULL when ks == n
    const int64_t *Uci_ks,       // the row index for nnz pattern of
                                 // U(k+1:ks-1,Q[ks])
    const int64_t *Ucx_ks,       // the value of i-th entry is found as
                                 // U->v[Uci[i]]->x[Ucx[i]]
    const int64_t Uc_ks_nz,      // # of nnz in U(k+1:ks-1,Q[ks])
    const int64_t k,             // current column index 0 <= k < n
    const int64_t ks,            // index of the diagonal to be swapped with,
                                 // k < ks <= n
    const SPEX_options *option   // command options
);

// perform diagonal permutation pivot update
SPEX_info spex_update_dppu1
(
    SPEX_matrix *L,              // matrix L
    SPEX_matrix *U,              // matrix U
    SPEX_matrix *rhos,           // array of scaled pivots
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,              // the index of first off-diag entry in
                                 // column k of L
    int64_t *h,                  // allocated vector that can be used for
                                 // history vector.  All entries are maintained
                                 // to be >= -1
    int64_t *Q,                  // column permutation
    int64_t *Q_inv,              // inverse of column permutation
    int64_t *P,                  // row permutation
    int64_t *P_inv,              // inverse of row permutation
    const int64_t k,             // current column index 0 <= k < n
    const int64_t ks,            // index of the diagonal to be swapped with,
                                 // k < ks <= n
    const SPEX_options *option   // command options
);

// perform diagonal permutation pivot update
SPEX_info spex_update_dppu2
(
    SPEX_matrix *L,              // matrix L
    SPEX_matrix *U,              // matrix U
    SPEX_matrix *rhos,           // array of scaled pivots
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *jnext,              // the index of first off-diag entry in
                                 // row k of U
    int64_t *h,                  // allocated vector that can be used for
                                 // history vector.  All entries are maintained
                                 // to be >= -1
    int64_t *Q,                  // column permutation
    int64_t *Q_inv,              // inverse of column permutation
    int64_t *P,                  // row permutation
    int64_t *P_inv,              // inverse of row permutation
    const int64_t k,             // current column index 0 <= k < n
    const int64_t ks,            // index of the diagonal to be swapped with,
                                 // k < ks <= n
    const SPEX_options *option   // command options
);

// perform history update for entries that would be in L and insert entries
// that would in U to corresponding row of U.
SPEX_info spex_update_finalize_and_insert_vk
(
    spex_scattered_vector *vk_dense, //scattered version of the solution for
                                 // LDx=v using the first k-1 columns of L
    int64_t *h,                  // history vector for vk_dense
    SPEX_matrix *U,              // matrix U
    SPEX_matrix *L,              // matrix L
    const SPEX_matrix *rhos,     // array of scaled pivots
    const int64_t *Q,            // the column permutation
    const int64_t *P_inv,        // inverse of row permutation
    const int64_t k,             // the column index in L that vk_dense
                                 // will be inserted
    const int64_t diag,          // the index of entry in vk_dense that
                                 // will be diagonal
    const SPEX_options *option   // command options
);

// insert an entry vi who has no pending scale to a scaled vector v, all
// v->x[i] will be scaled and S will be 1 after vi is inserted.
SPEX_info spex_update_insert_new_entry
(
    mpz_t vi,          // the entry to be inserted as i-th entry of v1
    SPEX_vector *v,    // the vector that would add new entry
    mpq_t S,           // pending scale for v1
    const int64_t i,   // the index of vi when inserted to v1
    const SPEX_options *option// command options
);

// perform one iteration of IPGE and perform skipped any scaling process.
SPEX_info spex_update_ipge // perform IPGE on x based on v
(
    spex_scattered_vector *sv_x,// array of size n for x in the scattered form.
                    // x could be dense by setting sv_x->i = NULL.
    int64_t *h,     // history vector for x, x[i] was last updated in the
                    // SPEX_FLIP(h[i])-th iteration
    int64_t *prev,  // prev is the index of the found previous entry of the last
                    // one (i.e., 2nd last entry) in v(perm). update if !prev
    SPEX_vector *v, // v is the vector that contains the j-th pivot
                    // used to compute x in the j-th IPGE iteration, which is
                    // the vector v in the equations mentioned above
    const int64_t *perm, // permutation
    const int64_t *perm_inv, // inverse of permutation
    const SPEX_matrix *rhos,// array of scaled pivots
    const int64_t j // column index of v
);

// perform REF triangular solve up to specified iteration, additional history
// update for certain entries should be done after calling this function.
SPEX_info spex_update_triangular_solve // perform REF triangular solve for LDx=v
(
    spex_scattered_vector *sv_x,// the scattered version of solution for LDx=v,
                        // using the first k-1 columns of L
    int64_t *x_top,     // P_inv[sv_x->i[0...(*x_top)]] <= (*last_update), that
                        // is, sv_x->i[0...(*x_top)] give the indices of all
                        // entries that are up-to-date. However, this is updated
                        // only when i_2ndlast is requested.
    int64_t *h,         // history vector for x
    int64_t *last_update,// the number of finished IPGE iterations, which is
                        // also the number of columns in L used last time
    int64_t *i_2ndlast, // i_2ndlast is the index of the found last nnz entry
                        // of x[P] less than n, this could be NULL if not needed
    const int64_t k,    // compute x up to k-th IPGE iteration, that is, using
                        // the first k-1 columns of L
    const SPEX_matrix *L,  // matrix L
    const SPEX_matrix *U,  // matrix U
    const SPEX_matrix *rhos,// array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
);

// sparse forward substitution, i.e., compute x = (LD)\v
SPEX_info spex_update_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    const SPEX_matrix *L,  // matrix L
    const int64_t *P,   // row permutation
    const SPEX_matrix *rhos,// array of scaled pivots
    int64_t *h          // history vector for x
);

// sparse REF backward substitution i.e., compute x = U\b
SPEX_info spex_update_backward_sub  // performs sparse REF backward substitution
(
    SPEX_vector *x,         // right hand side vector
    const SPEX_matrix *U,      // input upper triangular matrix
    const SPEX_matrix *rhos,// array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv    // inverse of column permutation
);


SPEX_info spex_update_solve_internal
(
    // Output
    SPEX_matrix **x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    const SPEX_matrix *b,   // a m*n dense matrix contains the right-hand-side
                            // vector
    SPEX_factorization *F,  // The SPEX LU or Cholesky factorization
    const bool transpose,   // whether computing Ax=b or ATx=b
    const SPEX_options* option // Command options
);
//------------------------------------------------------------------------------
// convert matrix between updatable and non-updatable format
//------------------------------------------------------------------------------

// NOTE: The L and U factorization from SPEX_Left_LU or the L from
// SPEX_Cholesky are all in SPEX_CSC format, and their columns and rows are
// permuted to be the same as the permuted matrix A(P,Q), and thus
// A(P,Q)=LD^(-1)U. However, all Update functions requires A=LD^(-1)U. This is
// the function to perform the coversion for L and U so that it meet the
// requirement for Update function and vice versa. Although the exact steps are
// slightly different, the overall result are equivalent to the following:
//
// To get updatable (dynamic_CSC MPZ) from non-updatable (CSC MPZ) matrix:
//     1. performing SPEX_matrix_copy to get the L and/or UT
//     SPEX_DYNAMIC_CSC form.
//     2. permute row indices of L or UT such that A=LD^(-1)U, or equivalently
//     A_out->v[j]->i[p] = perm[A_in->v[j]->i[p]], where A is either L or U.
//     3. canonicalize a SPEX_DYNAMIC_CSC matrix such that each column of the
//     input matrix have corresponding pivot as the first entry.
//
// To get non-updatable (CSC MPZ) from updatable (dynamic_CSC MPZ) matrix:
//     1. the inverse of the permutation
//     2. call SPEX_matrix_copy to obtain L and/or U in the SPEX_CSC format.

/*
SPEX_info spex_update_matrix_convert
(
    SPEX_factorization *F,       // converted CSC matrix
                                 // identity matrix if input as NULL
    const bool convertL,         // true if B->v[i] is the i-th col of B.
                                 // Otherwise, B->v[i] is the i-th row of B
    const SPEX_options *option
);*/
#endif

