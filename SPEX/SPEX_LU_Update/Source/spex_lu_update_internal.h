//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_lu_update_internal: include file for internal use in SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_LU_Update.h instead.

// uncommemd to help debugging
//#define SPEX_DEBUG


#ifndef SPEX_LU_UPDATE_INTERNAL_H
#define SPEX_LU_UPDATE_INTERNAL_H

#include "spex_util_internal.h"
#include "SPEX_LU_Update.h"
#ifdef SPEX_LU_UPDATE_TCOV
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

/*
SPEX_info spex_cast_array
(
    void *Y,                // output array, of size n
    SPEX_type ytype,        // type of Y
    void *X,                // input array, of size n
    SPEX_type xtype,        // type of X
    int64_t n,              // size of Y and X
    mpq_t y_scale,          // scale factor applied if y is mpz_t
    mpq_t x_scale,          // scale factor applied if x is mpz_t
    const SPEX_options *option
) ;

SPEX_info spex_cast_matrix
(
    SPEX_mat **Y_handle,     // nz-by-1 dense matrix to create
    SPEX_type Y_type,           // type of Y
    SPEX_mat *A,             // matrix with nz entries
    const SPEX_options *option
) ;

// (void *) pointer to the values of A.  A must be non-NULL with a valid type
#define SPEX_X(A)                                                           \
    ((A->type == SPEX_MPZ  ) ? (void *) A->x.mpz   :                        \
    ((A->type == SPEX_MPQ  ) ? (void *) A->x.mpq   :                        \
    ((A->type == SPEX_MPFR ) ? (void *) A->x.mpfr  :                        \
    ((A->type == SPEX_INT64) ? (void *) A->x.int64 : (void *) A->x.fp64))))


// return an error if A->kind (csc, triplet, dense) is wrong
#define SPEX_REQUIRE_KIND(A,required_kind) \
    if (A == NULL || A->kind != required_kind) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_KIND(A,required_kind) \
    ASSERT (A != NULL && A->kind == required_kind)

// return an error if A->type (mpz, mpq, mpfr, int64, or double) is wrong
#define SPEX_REQUIRE_TYPE(A,required_type) \
    if (A == NULL || A->type != required_type) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_TYPE(A,required_type) \
    ASSERT (A != NULL && A->type == required_type)

// return an error if A->kind or A->type is wrong
#define SPEX_REQUIRE(A,required_kind,required_type)     \
    SPEX_REQUIRE_KIND (A,required_kind) ;               \
    SPEX_REQUIRE_TYPE (A,required_type) ;

#define ASSERT_MATRIX(A,required_kind,required_type)    \
    ASSERT_KIND (A,required_kind) ;                     \
    ASSERT_TYPE (A,required_type) ;
*/

/*
typedef struct
{
    int64_t n;    // number of entries
    int64_t nz;   // number of nonzeros
    int64_t *i;   // array of size nzmax that contains the column/row indices
                  // of nnz
    //bool i_shallow;// if true, i is shallow, and original i should be updated
                  // when this struct is free'd
    mpz_t *x;     // array of size n that contains the values of each entry
} spex_scattered_vector;

spex_scattered_vector* spex_create_scattered_vector
(
    const int64_t n             // number of entries in v
);

void spex_delete_scattered_vector
(
    spex_scattered_vector **sv  // scattered vector to be deleted
);*/
#define spex_scattered_vector SPEX_vector
#define spex_scattered_vector_alloc SPEX_vector_alloc
#define spex_scattered_vector_free SPEX_vector_free

SPEX_info spex_get_scattered_v
(
    spex_scattered_vector **sv_handle,// output vector in scattered form
    SPEX_vector *v,              // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const bool keep_v            // indicate if the mpz values should be kept
);

SPEX_info spex_cppu
(
    SPEX_mat *L,     // matrix L
    SPEX_mat *U,     // matrix U
    SPEX_matrix *S,  // a 2*n dense mpq matrix that stores pending scales
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                     // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k, // current column index 0 <= k < n
    const int64_t ks // index of the diagonal to be swapped with, [0,n)
);

SPEX_info spex_dppu1
(
    SPEX_mat *L,     // matrix L
    SPEX_mat *U,     // matrix U
    SPEX_matrix *S,  // a 2*n dense mpq matrix that stores pending scales
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


SPEX_info spex_dppu2
(
    SPEX_mat *L,     // matrix L
    SPEX_mat *U,     // matrix U
    SPEX_matrix *S,  // a 2*n dense mpq matrix that stores pending scales
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

SPEX_info spex_get_nnz_pattern    // find the nnz pattern of L and U
(
    // OUTPUT:
    int64_t **Lr_offdiag,         // Lr_offdiag[k] gives the column index of the
                                  // last off-diagonal nnz in k-th row of L.
                                  // -1 if no off diagonal entry
    int64_t **Uc_offdiag,         // Lr_offdiag[k] gives the column index of the
    int64_t **Uci,                // the row index for col-wise nnz pattern of U
    int64_t **Ucp,                // col pointers for col-wise pattern of U
    int64_t **Ucx,                // find the value of k-th entry as
                                  // U->v[Uci[k]]->x[Ucx[k]]
    // INPUT:
    const SPEX_mat *L,            // the target matrix L
    const SPEX_mat *U,            // the target matrix U
    const int64_t *P_inv,         // inverse of row permutation
    const int64_t *Q_inv,         // inverse of col permutation
    const int64_t k,              // index of column to be replaced
    const SPEX_options *option    // command option
);

SPEX_info spex_finalize_and_insert_vk
(
    spex_scattered_vector *vk_dense, //scattered version of the solution for
                      // LDx=v using the first k-1 columns of L
    int64_t *h,       // history vector for vk_dense
    SPEX_mat *U,      // matrix U
    SPEX_mat *L,      // matrix L
    SPEX_matrix *S,   // a 2*n dense mpq matrix that stores pending scales
    const mpz_t *sd,  // array of scaled pivots
    const int64_t *Q, // the column permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t k,  // the column index in L that vk_dense will be inserted
    const int64_t diag,// the index of entry in vk_dense that will be diagonal
    const mpq_t one   // a const mpq number, just to avoid constantly alloc
);

SPEX_info spex_find_next_nz
(
    int64_t *next,                  // the col/row index of next nnz
    spex_scattered_vector *Ak_dense,// the scattered vector
    int64_t *perm_inv,              // inverse of permutation
    int64_t k
);

SPEX_info spex_insert_new_entry
(
    mpz_t vi,          // the entry to be inserted as i-th entry of v1
    SPEX_vector *v,    // the vector that would add new entry
    mpq_t S,           // pending scale for v1
    const int64_t i,   // the index of vi when inserted to v1
    const mpq_t one    // a constant mpq number, just to avoid constantly alloc
);

SPEX_info spex_ipge // perform IPGE on x based on v
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
    mpq_t v_scale1, // the first pending scale for v
    const int64_t j // column index of v
);

SPEX_info spex_triangular_solve // perform REF triangular solve for LDx=v
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
    SPEX_matrix *S,     // a 2*n dense mpq matrix that stores pending scales
    const mpz_t *sd,    // array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
);

SPEX_info spex_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    const SPEX_mat *L,  // matrix L
    const int64_t *P,   // row permutation
    const mpz_t *sd,    // array of scaled pivots
    SPEX_matrix *S,     // a 2*n dense mpq matrix that stores pending scales
    int64_t *h,          // history vector for x
    const bool Is_trans    // true if solving A'x = b
);

SPEX_info spex_backward_sub  // performs sparse REF backward substitution
(
    SPEX_vector *x,         // right hand side vector
    const SPEX_mat *U,      // input upper triangular matrix
    const SPEX_matrix *S,   // a 2*n dense mpq matrix that stores pending scales
    const mpz_t *sd,        // array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv,    // inverse of column permutation
    const bool Is_trans    // true if solving A'x = b
);

SPEX_info spex_verify
(
    bool *correct,         // indicate if the verification is passed
    const SPEX_mat *L,     // lower triangular matrix
    const SPEX_mat *U,     // upper triangular matrix
    const SPEX_mat *A,     // Input matrix
    int64_t *h,            // history vector
    SPEX_matrix *S,        // a 2*n dense mpq matrix that stores pending scales
    const mpz_t *sd,       // array of scaled pivots
    const int64_t *P,      // row permutation
    const int64_t *Q_inv,  // inverse of column permutation
    const SPEX_options *option// command options
);

#endif

