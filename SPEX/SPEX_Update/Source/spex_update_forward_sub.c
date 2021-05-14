//------------------------------------------------------------------------------
// SPEX_Update/spex_update_forward_sub.c: sparse forward substitution,
// i.e., compute x = (LD)\v
//------------------------------------------------------------------------------
    
// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------
    
// Purpose: This function is to perform sparse forward substitution, which is
// essentially the same as sparse REF triangular solve for LDx=v, but with v
// as a dense vector. This function assumes v is in the same row permutation as
// L. This function take v as input using x_input, and the solution is stored
// in x_output. In case when v has multiple column, simply iterate this function
// and solve each column at each iteration.

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)

SPEX_info spex_update_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    const SPEX_matrix *L,  // matrix L
    const int64_t *P,   // row permutation
    const SPEX_matrix *rhos,// array of scaled pivots
    int64_t *h          // history vector for x
)
{
    // TODO only perform input check in user callable functions
    if (!x || !h || !L ||!P || !rhos)
    {
        return SPEX_INCORRECT_INPUT;
    }

    SPEX_info info;
    int sgn;
    int64_t i, n = L->n;
    mpq_t x_scale; SPEX_MPQ_SET_NULL(x_scale);// TODO delete
    SPEX_CHECK(SPEX_mpq_init(x_scale));
    SPEX_CHECK(SPEX_mpq_set_ui(x_scale, 1, 1));
    for (i = 0; i < n; i++)
    {
        h[i] = -1; 
    }

    for (i = 0; i < n; i++)
    {
        // skip if x(P[i]) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x[P[i]]));
        if (sgn == 0)       { continue; }
        //printf("%ld-th IPGE\n",i);

        // check if the first entry of L->v[i] is the corresponding pivot
        ASSERT(L->v[i]->i[0] == P[i]);
        // perform i-th IPGE update for x
#if 1
        SPEX_CHECK(spex_update_ipge(x, h, NULL, L->v[i], P, NULL, rhos, i));
#else
        /*
        SPEX_CHECK(spex_ipge1(x, h, NULL, L, P, NULL, rhos, SL(i), i));
        */
        
        /*
        // incorrect
        SPEX_CHECK(spex_ipge2(x, h, NULL, L, P, NULL, i));
        if (i-1 > 0)
        {
            SPEX_CHECK(SPEX_mpz_mul(x->x[P[i]],
                                    x->x[P[i]], SPEX_MPQ_NUM(SL(i-1))));
            SPEX_CHECK(SPEX_mpz_divexact(x->x[P[i]],
                                    x->x[P[i]], SPEX_MPQ_DEN(SL(i-1))));
        }
        */
        SPEX_CHECK(spex_ipge3(x, h, NULL, L, P, NULL, rhos, x_scale,SL(i), i));
        //gmp_printf("%Qd\n",x_scale);
#endif
    }

    return SPEX_OK; 
}
