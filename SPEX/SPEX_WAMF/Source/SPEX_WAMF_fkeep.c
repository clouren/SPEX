//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_fkeep: Drop entries for which fkeep(A(i,j)) is false
//------------------------------------------------------------------------------
// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 
 */

int64_t SPEX_WAMF_fkeep 
(
    SPEX_matrix *A, 
    int64_t (*fkeep) (int64_t, int64_t, double, void *), 
    void *other
)
{
    int64_t j, p, nz = 0, n, *Ap, *Ai ;
    double *Ax ;
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x.fp64 ;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
            if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
            {
                if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Ai [nz++] = Ai [p] ;
            }
        }
    }
    Ap [n] = nz ;                           /* finalize A */
    return (nz) ;
}
