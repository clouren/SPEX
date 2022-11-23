//------------------------------------------------------------------------------
// SPEX_Update/spex_update_insert_new_entry: insert an entry vi who has no
// pending scale to a scaled vector v, all v->x[i] will be scaled and S will be
// 1 after vi is inserted.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is call to insert vi as i-th entry of vector v, where
// v could be column of L or row of U in frame j. S is the pending scales for
// vector v. Since vi has no pending scale factor and v has pending scale
// factor, we will apply S to all v->x[i] before insert vi into v. S is reset to
// 1 when finished.
//
// NOTE: the value of vi is inserted by swapping with corresponding mpz_t value
// in v, which should be considered as undefined after this function.

#include "spex_update_internal.h"

SPEX_info spex_update_insert_new_entry
(
    mpz_t vi,           // the entry to be inserted as i-th entry of v
    SPEX_vector v,      // the vector that would add new entry
    mpq_t S,            // pending scale for v
    const int64_t i,    // the index of vi when inserted to v
    const SPEX_options option
)
{

    SPEX_info info;
    int r;
    int64_t p;

    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, S, 1, 1));
    if (r != 0) //S != 1
    {
        for (p = 0; p < v->nz; p++)
        {
            // Since entries in v will be integer after scale, we can
            // perform division first to make it small, and this
            // division will preserve integer propety
            SPEX_CHECK(SPEX_mpz_divexact(v->x[p], v->x[p], SPEX_MPQ_DEN(S)));
            SPEX_CHECK(SPEX_mpz_mul(v->x[p], v->x[p], SPEX_MPQ_NUM(S)));
        }
        SPEX_CHECK(SPEX_mpq_set_ui(S, 1, 1));
    }
    // append vi to v
    if (v->nz == v->nzmax)
    {
        // reallocate the nonzero pattern if needed
        SPEX_CHECK(SPEX_vector_realloc(v, v->nzmax+1, option));
    }
    v->i[v->nz] = i;
    SPEX_CHECK(SPEX_mpz_swap(v->x[v->nz], vi));
    v->nz ++;

    return SPEX_OK;
}
