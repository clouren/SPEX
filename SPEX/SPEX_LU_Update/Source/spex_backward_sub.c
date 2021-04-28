//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_backward_sub: sparse REF backward substitution (x = U\b)
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF backward substitution, solving
 * the system Ux = b. x is internally multiplied by the determinant of A to
 * maintain integral. In addition, x is also permuted by Q. Therefore the
 * real solution should be found as x(i) = x(Q_inv(i))/det(A).
 *
 * U is a sparse mpz matrix, which is stored by row, and x is a dense mpz
 * vector.
 *
 * The input argument x contains b on input, and it is overwritten on output
 * by the solution x.
 */

#define SPEX_FREE_ALL                \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_lu_update_internal.h"

#define SL(k) S->x.mpq[2*(k)]
#define SU(k) S->x.mpq[1+2*(k)]

SPEX_info spex_backward_sub  // performs sparse REF backward substitution
(
    SPEX_vector *x,         // right hand side vector
    const SPEX_mat *U,      // input upper triangular matrix
    const SPEX_matrix *S,   // a 2*n dense mpq matrix that stores pending scales
    const mpz_t *sd,        // array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv,    // inverse of column permutation
    const bool Is_trans    // true if solving A'x = b
)
{
    SPEX_info info ;
    int sgn;
    int64_t i, j, real_i, real_j, p, n = U->n;
    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    // Start at x[n-1], since x[n] will remain the same
    for (i = n-2; i >= 0; i--)
    {
        // tmpz = 0
        SPEX_CHECK(SPEX_mpz_set_ui(tmpz, 0));
        real_i = P[i];

        // skip the diagonal (pivot entry), which locates at p = 0
        for (p = 1; p < U->v[i]->nz; p++) // i-th row of U
        {
            j = U->v[i]->i[p];// the real col index is Q_inv[j]
            real_j = P[Q_inv[j]];

            // skip if corresponding entry in x is zero
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x[real_j]));
            if (sgn == 0) {continue;}

            // tmpz -= U(i,j)*x[j]
            SPEX_CHECK(SPEX_mpz_submul(tmpz, U->v[i]->x[p],x->x[real_j]));
        }
        // tmpz = tmpz*S(2,i)
        if (Is_trans)
        {
            SPEX_CHECK(SPEX_mpz_divexact(tmpz, tmpz, SPEX_MPQ_DEN(SL(i))));
            SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(SL(i))));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_divexact(tmpz, tmpz, SPEX_MPQ_DEN(SU(i))));
            SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(SU(i))));
        }

        // x[i] = x[i]*sd[n-1]+tmpz
        SPEX_CHECK(SPEX_mpz_mul(x->x[real_i], x->x[real_i], sd[n-1]));
        SPEX_CHECK(SPEX_mpz_add(x->x[real_i], x->x[real_i], tmpz));

        // x[i] = x[i]/sd[i]
        SPEX_CHECK(SPEX_mpz_divexact(x->x[real_i], x->x[real_i], sd[i]));
    }

    SPEX_FREE_ALL;
    return (SPEX_OK) ;
}
