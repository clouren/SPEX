//------------------------------------------------------------------------------
// SPEX_Update/spex_update_insert_new_entry.c: insert an entry vi who has no
// pending scale to a scaled vector v, find a common scale factor for v and vi
// and update (scale up) entries in v properly.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is call to insert vi as i-th entry of vector v, where
// v could be column of L or row of U in frame j. S is the pending scales for
// vector v. Since vi has no pending scale factor and v has pending scale
// factor, we will find a new common factor for vi and v, and properly scale v
// and vi before insert vi into v.
//
// NOTE: the value of vi is inserted by swapping with corresponding mpz_t value
// in v, which should be considered as undefined after this function.

#include "spex_update_internal.h"

SPEX_info spex_update_insert_new_entry
(
    mpz_t vi,          // the entry to be inserted as i-th entry of v
    SPEX_vector *v,   // the vector that would add new entry
    mpq_t S,          // pending scale for v
    const int64_t i,   // the index of vi when inserted to v
    const mpq_t one    // a constant mpq number, just to avoid constantly alloc
)
{
    SPEX_info info;
    int r;
    int64_t p;

    SPEX_CHECK(SPEX_mpq_equal(&r, S, one));
    if (r == 0) //S != 1
    {
#if 0
        // find the gcd of inserted entry and numerator of S
        SPEX_CHECK(SPEX_mpz_gcd(gcd, vi, SPEX_MPQ_NUM(S)));
        SPEX_CHECK(SPEX_mpz_divexact(vi, vi, gcd));
        SPEX_CHECK(SPEX_mpz_divexact(SPEX_MPQ_NUM(S), SPEX_MPQ_NUM(S), gcd));
#endif
        for (p = 0; p < v->nz; p++)
        {
            // Since entries in v will be integer after scale, we can
            // perform division first to make it small, and this
            // division will preserve integer propety
            SPEX_CHECK(SPEX_mpz_divexact(v->x[p], v->x[p], SPEX_MPQ_DEN(S)));
            SPEX_CHECK(SPEX_mpz_mul(v->x[p], v->x[p], SPEX_MPQ_NUM(S)));
        }
#if 0
        SPEX_CHECK(SPEX_mpq_set_z(S, gcd));
#else
        SPEX_CHECK(SPEX_mpq_set(S, one));
#endif
    }
    // append vi to v
    if (v->nz == v->nzmax)
    {
        // reallocate the nonzero pattern if needed
        SPEX_CHECK(SPEX_vector_realloc(v, v->nzmax+1));
    }
    v->i[v->nz] = i;
    SPEX_CHECK(SPEX_mpz_swap(v->x[v->nz], vi));
    v->nz ++;

    return SPEX_OK;
}
