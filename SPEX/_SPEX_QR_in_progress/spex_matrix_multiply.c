//------------------------------------------------------------------------------
// SPEX_Utilities/spex_matrix_multiply: Multiply two matrices
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020-2023, Lorena Mejia Domenzain, Christopher Lourenco,
// Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

//TODO
SPEX_info spex_matrix_multiply 
(
    //Output
    SPEX_matrix C,
    //Input
    const SPEX_matrix A, 
    const SPEX_matrix B
)
{
    //Declare variables
    csi p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    double *x, *Bx, *Cx ;
    //Allocate memory for C
    cs *C ;

    //Check inputs
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;    
    if (A->n != B->m) return (NULL) ;


    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; bnz = B->p [n] ;
    w = cs_calloc (m, sizeof (csi)) ;                    /* get workspace */
    values = (A->x.mpz != NULL) && (B->x.mpz != NULL) ;
    x = values ? cs_malloc (m, sizeof (double)) : NULL ; /* get workspace */
    C = cs_spalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (cs_done (C, w, x, 0)) ;             /* out of memory */
        } 
        C->p [j] = nz ;                   /* column j of C starts here */
        for (p = B->p [j] ; p < B->p [j+1] ; p++)
        {
            nz = spex_scatter (A, B->i [p], B->x ? B->x [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = C->p [j] ; p < nz ; p++) C->x [p] = x [C->i [p]] ;
    }
    C->p [n] = nz ;                       /* finalize the last column of C */

    FREE_WORKSPACE;
    return SPEX_OK;
}