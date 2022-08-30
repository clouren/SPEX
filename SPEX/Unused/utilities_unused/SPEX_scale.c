//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_scale: FIXME
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose:  FIXME
 *
 * Input/Output arguments:
 *
 */

#define SPEX_FREE_WORKSPACE                 \
{                                            \
    SPEX_MPQ_CLEAR(scale);                   \
}

#include "spex_util_internal.h"

SPEX_info SPEX_scale
(
    // Output
    SPEX_matrix* x,
    // Input
    const mpq_t scaling_num, //numerator
    const mpq_t scaling_den, //denominator
    const SPEX_options* option        // command options
)
{
    SPEX_info info;

    int i, nz;

    // Scale is the scaling factor for the solution vectors.
    // When the forward/backsolve is complete, the entries in
    // x are rational, but are solving the scaled linear system
    // A' x = b' (that is if A had input which was rational or floating point
    // and had to be converted to integers). Thus, the scale here is used
    // to convert x into into the actual solution of A x = b. 
    mpq_t scale;
    SPEX_MPQ_SET_NULL(scale);

    SPEX_CHECK(SPEX_mpq_init(scale));

    // set the scaling factor scale = scaling_num / scaling_den
    SPEX_CHECK(SPEX_mpq_div(scale, scaling_num, scaling_den));

    // Apply scaling factor, but ONLY if it is different from a scaling factor
    // to 1
    int r;
    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, scale, 1, 1));
    if (r != 0)
    {
        nz = x->m * x->n;
        for (i = 0; i < nz; i++)
        {
            SPEX_CHECK(SPEX_mpq_mul(x->x.mpq[i], x->x.mpq[i], scale));
        }
    }

  //--------------------------------------------------------------------------
  // Free all workspace and return success
  //--------------------------------------------------------------------------
  SPEX_FREE_WORKSPACE;
  return SPEX_OK;
}
