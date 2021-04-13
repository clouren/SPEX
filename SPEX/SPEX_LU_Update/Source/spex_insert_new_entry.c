//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_insert_new_entry.c: insert an entry vi who has no pending
// scale to a scaled vector v1, find a common scale factor for v1 and vi and
// update (scale up) entries in v1 properly.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is call to insert vi as i-th entry of vector v1, which
// could be column of L or row of U in frame j. v1 and v2 together constitute
// frame j of the frame matrix. S1, S2 and S3 are correspondingly the pending
// scales for vector v1, vector v2, and all entries in frame j (If v1 is the
// j-th column of L, then S1, S2, S3 are scales for j-th column of L, j-th row
// of U and j-th frame of the frame matrix). Since vi has no pending scale
// factor and v1 has pending scale factor, we will find a new common factor for
// vi and v1, and properly scale v1 and vi before insert vi into v1.
//
// NOTE: the value of vi is inserted by swapping with corresponding mpz_t value
// in v1, which should be considered as undefined after this function.

#include "spex_lu_update_internal.h"

SPEX_info spex_insert_new_entry
(
    mpz_t vi,          // the entry to be inserted as i-th entry of v1
    SPEX_vector *v1,   // the vector that would add new entry
    mpq_t S1,          // pending scale for v1
    const SPEX_vector *v2,// the other vector that is in same frame as v1
    mpq_t S2,          // pending scale for v2
    mpq_t S3,          // pending scale for frame that holds v1 and v2
    mpz_t d,           // the unscale pivot in frame of v1 and v2
    const int64_t i,   // the index of vi when inserted to v1
    const int64_t v2_diag,// the pointer to the diagonal entry in v2
    const mpq_t one    // a constant mpq number, just to avoid constantly alloc
)
{
    SPEX_info info;
    int r;
    int64_t p;

    // d = v2(v2_diag) = d*S1
    SPEX_CHECK(SPEX_mpz_set(d, v2->x[v2_diag]));
    // mpq_equal is said to be faster than mpq_cmq
    SPEX_CHECK(SPEX_mpq_equal(&r, one, one));
    SPEX_CHECK(SPEX_mpq_equal(&r, S3, one));
    if (r == 0) // S3 != 1
    {
#if 0
        // find the gcd of inserted entry vi and numerator of scaling
        // factor for frame i
        SPEX_CHECK(SPEX_mpz_gcd(gcd, vi, SPEX_MPQ_NUM(S3)));
        SPEX_CHECK(SPEX_mpz_divexact(vi, vi, gcd));
        SPEX_CHECK(SPEX_mpz_divexact(SPEX_MPQ_NUM(S3), SPEX_MPQ_NUM(S3), gcd));
#endif
        // [S1;S2;S3] = [S1*S3; S2*S3; 1]
        SPEX_CHECK(SPEX_mpq_mul(S1, S1, S3));
        SPEX_CHECK(SPEX_mpq_mul(S2, S2, S3));
#if 0
        SPEX_CHECK(SPEX_mpq_set_z(S3, gcd));
#else
        SPEX_CHECK(SPEX_mpq_set(S3, one));
#endif
    }

    SPEX_CHECK(SPEX_mpq_equal(&r, S1, one));
    if (r == 0) //S1 != 1
    {
#if 0
        // find the gcd of inserted entry and numerator of S1
        SPEX_CHECK(SPEX_mpz_gcd(gcd, vi, SPEX_MPQ_NUM(S1)));
        SPEX_CHECK(SPEX_mpz_divexact(vi, vi, gcd));
        SPEX_CHECK(SPEX_mpz_divexact(SPEX_MPQ_NUM(S1), SPEX_MPQ_NUM(S1), gcd));
#endif
        for (p = 0; p < v1->nz; p++)
        {
            // Since entries in v1 will be integer after scale, we can
            // perform division first to make it small, and this
            // division will preserve integer propety
            SPEX_CHECK(SPEX_mpz_divexact(v1->x[p], v1->x[p],
                                       SPEX_MPQ_DEN(S1)));
            SPEX_CHECK(SPEX_mpz_mul(v1->x[p], v1->x[p], SPEX_MPQ_NUM(S1)));
        }
#if 0
        SPEX_CHECK(SPEX_mpq_set_z(S1, gcd));
#else
        SPEX_CHECK(SPEX_mpq_set(S1, one));
#endif
    }
    // append vi to v1
    if (v1->nz == v1->nzmax)
    {
        // reallocate the nonzero pattern if needed
        SPEX_CHECK(SPEX_vector_realloc(v1, 2*(v1->nzmax)));
    }
    v1->i[v1->nz] = i;
    SPEX_CHECK(SPEX_mpz_swap(v1->x[v1->nz], vi));
    v1->nz ++;

    return SPEX_OK;
}
