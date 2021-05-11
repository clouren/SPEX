//------------------------------------------------------------------------------
// SPEX_Update/spex_update_internal: include file for internal use in SPEX_Update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_Update.h instead.

// uncommemd to help debugging
#define SPEX_DEBUG


#ifndef SPEX_UPDATE_INTERNAL_H
#define SPEX_UPDATE_INTERNAL_H

#include "spex_util_internal.h"
#include "SPEX_Update.h"
#ifdef SPEX_UPDATE_TCOV
    /* include this header to make macro GOTCHA available */
    #include "../Tcov/tcov_malloc_test.h"
#endif
#ifndef GOTCHA
    #define GOTCHA
#endif

#ifdef SPEX_DEBUG
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
#endif


// ============================================================================
//                           Internal Functions
// ============================================================================

// check if SPEX_initialize* has been called
/*bool spex_initialized ( void ) ;        // true if called, false if not
void spex_set_initialized (bool s) ;    // set global initialzed flag to s
*/

#define spex_scattered_vector SPEX_vector
#define spex_scattered_vector_alloc SPEX_vector_alloc
#define spex_scattered_vector_free SPEX_vector_free

SPEX_info spex_update_get_scattered_v
(
    spex_scattered_vector **sv_handle,// output vector in scattered form
    SPEX_vector *v,              // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const bool keep_v            // indicate if the mpz values should be kept
);

SPEX_info spex_update_cppu
(
    SPEX_mat *L,     // matrix L
    SPEX_mat *U,     // matrix U
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation (unchanged on output)
    int64_t *P_inv,  // inverse of row permutation (unchanged on output)

    // the col-wise nnz pattern of U can be NULL when ks == n
    const int64_t *Uci_ks,// the row index for nnz pattern of U(k+1:ks-1,Q[ks])
    const int64_t *Ucx_ks,// the value of i-th entry is found as
                     // U->v[Uci[i]]->x[Ucx[i]]
    const int64_t Uc_ks_nz,// # of nnz in U(k+1:ks-1,Q[ks])
    const int64_t k, // current column index 0 <= k < n
    const int64_t ks // index of the diagonal to be swapped with, [0,n)
);

SPEX_info spex_update_dppu1
(
    SPEX_mat *L,     // matrix L
    SPEX_mat *U,     // matrix U
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
);


SPEX_info spex_update_dppu2
(
    SPEX_mat *L,     // matrix L
    SPEX_mat *U,     // matrix U
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
);

SPEX_info spex_update_finalize_and_insert_vk
(
    spex_scattered_vector *vk_dense, //scattered version of the solution for
                      // LDx=v using the first k-1 columns of L
    int64_t *h,       // history vector for vk_dense
    SPEX_mat *U,      // matrix U
    SPEX_mat *L,      // matrix L
    const mpz_t *sd,  // array of scaled pivots
    const int64_t *Q, // the column permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t k,  // the column index in L that vk_dense will be inserted
    const int64_t diag,// the index of entry in vk_dense that will be diagonal
    const mpq_t one   // a const mpq number, just to avoid constantly alloc
);

SPEX_info spex_update_find_next_nz
(
    int64_t *next,                  // the col/row index of next nnz
    spex_scattered_vector *Ak_dense,// the scattered vector
    int64_t *perm_inv,              // inverse of permutation
    int64_t k
);

SPEX_info spex_update_insert_new_entry
(
    mpz_t vi,          // the entry to be inserted as i-th entry of v1
    SPEX_vector *v,    // the vector that would add new entry
    mpq_t S,           // pending scale for v1
    const int64_t i,   // the index of vi when inserted to v1
    const mpq_t one    // a constant mpq number, just to avoid constantly alloc
);

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
    const mpz_t *sd,// array of scaled pivots
    const int64_t j // column index of v
);

SPEX_info spex_update_triangular_solve // perform REF triangular solve for LDx=v
(
    spex_scattered_vector *sv_x,// the scattered version of solution for LDx=v,
                        // using the first k-1 columns of L
    int64_t *h,         // history vector for x
    int64_t *last_update,// the number of finished IPGE iterations, which is
                        // also the number of columns in L used last time
    int64_t *i_2ndlast, // i_2ndlast is the index of the found last nnz entry
                        // of x[P] less than n, this could be NULL if not needed
    const int64_t k,    // compute x up to k-th IPGE iteration, that is, using
                        // the first k-1 columns of L
    const SPEX_mat *L,  // matrix L
    const SPEX_mat *U,  // matrix U
    const mpz_t *sd,    // array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
);

SPEX_info spex_update_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    const SPEX_mat *L,  // matrix L
    const int64_t *P,   // row permutation
    const mpz_t *sd,    // array of scaled pivots
    int64_t *h          // history vector for x
);

SPEX_info spex_update_backward_sub  // performs sparse REF backward substitution
(
    SPEX_vector *x,         // right hand side vector
    const SPEX_mat *U,      // input upper triangular matrix
    const mpz_t *sd,        // array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv    // inverse of column permutation
);

SPEX_info spex_update_verify
(
    bool *correct,         // indicate if the verification is passed
    const SPEX_mat *L,     // lower triangular matrix
    const SPEX_mat *U,     // upper triangular matrix
    const SPEX_mat *A,     // Input matrix
    int64_t *h,            // history vector
    const mpz_t *sd,       // array of scaled pivots
    const int64_t *P,      // row permutation
    const int64_t *Q_inv,  // inverse of column permutation
    const SPEX_options *option// command options
);

#endif

