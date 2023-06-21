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
    int64_t p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    SPEX_matrix *x, *Bx, *Cx ;
    //Allocate memory for C
    SPEX_matrix *C;

    //Check inputs
    //if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;    
    if (A->n != B->m) return (NULL) ;


    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; bnz = B->p [n] ;
    w = cs_calloc (m, sizeof (int64_t)) ;                    /* get workspace */
    values = (A->x.mpz != NULL) && (B->x.mpz != NULL) ;
    //x = values ? cs_malloc (m, sizeof (double)) : NULL ; /* get workspace */
    if(values)
    {
        SPEX_CHECK(SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPZ, m, 1, m,
        false, true, NULL));
    }
    else
    {
        return SPEX_OUT_OF_MEMORY;
    }
    
    SPEX_CHECK(SPEX_matrix_allocate(&C, SPEX_CSC, SPEX_MPZ, m, n, anz+bnz,
        false, true, NULL));
    if (!C || !w || (values && !x)) return (SPEX_OUT_OF_MEMORY) ;
    
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (SPEX_OUT_OF_MEMORY) ;             /* out of memory */
        } 
        C->p [j] = nz ;                   /* column j of C starts here */
        for (p = B->p [j] ; p < B->p [j+1] ; p++)
        {
            spex_scatter (A, B->i [p], B->x.mpz ? B->x.mpz [p] : 1, w, x, j+1, C, nz) ;//returns nz
        }
        if (values) for (p = C->p [j] ; p < nz ; p++) C->x.mpz [p] = x->x.mpz [C->i [p]] ;
    }
    C->p [n] = nz ;                       /* finalize the last column of C */

    FREE_WORKSPACE;
    return SPEX_OK;
}
