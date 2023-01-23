//------------------------------------------------------------------------------
// SPEX_tripread: reads a matrix stored in triplet format
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function reads a double or mpz matrix stored in triplet format
 * This format used can be seen in any of the example mat files.
 *
 * The first line of the file contains three integers: m, n, nnz,
 * where the matrix is m-by-n with nnz entries.
 *
 * This is followed by nnz lines, each containing a single triplet: i, j, aij,
 * which defines the row index (i), column index (j), and value (aij) of
 * the entry A(i,j).  
 */

#include "demos.h"

SPEX_info SPEX_tripread
(
    SPEX_matrix *A_handle,      // Matrix to be populated
    FILE *file,                 // file to read from (must already be open)
    SPEX_type C_type,          // C->type: mpz_t, mpq_t, mpfr_t, int64_t, or double
    SPEX_options option
)
{
    SPEX_info info ;
    if (A_handle == NULL || file == NULL)
    {
        printf ("invalid input\n") ;
        return SPEX_INCORRECT_INPUT;
    }
    
    switch (C_type)
    {
        case SPEX_MPZ:
            SPEX_CHECK(SPEX_tripread_mpz(&A_handle, file, option));
            break;
        case SPEX_FP64:
            SPEX_CHECK(SPEX_tripread_double(&A_handle, file, option));
            break;
    }

    return SPEX_OK;
}