//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_forward_sub.c: sparse forward substitution (x = (LD)\v)
//------------------------------------------------------------------------------
    
// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------
    
// Purpose: This function is to perform sparse forward substitution, which is
// essentially the same as sparse REF triangular solve for LDx=v, but with v
// as a dense vector. This function assumes v is in the same row permutation as
// L. This function take v as input using x_input, and the solution is stored
// in x_output. In case when v has multiple column, simply iterate this function
// and solve each column at each iteration.

#include "spex_lu_update_internal.h"

#define SL(k) S->x.mpq[2*(k)]
#define SU(k) S->x.mpq[1+2*(k)]

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
)
{
    if (!x || !h || !L || !S || !P || !sd)
    {
        return SPEX_INCORRECT_INPUT;
    }

    SPEX_info info;
    int sgn;
    int64_t i, n = L->n;
    for (i = 0; i < n; i++)
    {
        h[i] = -1; 
    }

    for (i = 0; i < n; i++)
    {
        // skip if x(P[i]) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x[P[i]]));
        if (sgn == 0)       { continue; }

        // check if the first entry of L->v[i] is the corresponding pivot
        ASSERT(L->v[i]->i[0] == P[i]);
        // perform i-th IPGE update for x
        if (Is_trans)
        {
            SPEX_CHECK(spex_ipge(x, h, NULL, L->v[i], P, NULL, sd, SU(i), i));
        }
        else
        {
            SPEX_CHECK(spex_ipge(x, h, NULL, L->v[i], P, NULL, sd, SL(i), i));
        }
    }

    return SPEX_OK; 
}
