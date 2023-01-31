//------------------------------------------------------------------------------
// Demo/spex_demo_tripread: reads a matrix stored in triplet format of a given type
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function reads a matrix stored in triplet format of a given type
 * This format used is illustrated in the example mat files.
 *
 * The first line of the file contains three integers: m, n, nnz,
 * where the matrix is m-by-n with nnz entries.
 *
 * Each of the following nnz lines contains a single triplet: i, j, aij,
 * which defines the row index (i), column index (j), and value (aij) of
 * the entry A(i,j).  
 */

// Lore FIXME: We had a miscommunication error. What we need is that the other two
// trireads---which overlap about 80%, including all of the flow, logic and error 
// checking---should be merged into a single file. This eliminates the unecessary 
// and difficult to maintain redundancies.


#include "demos.h"

SPEX_info spex_demo_tripread
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
        printf ("invalid input\n");
        return SPEX_INCORRECT_INPUT;
    }
    
    switch (C_type)
    {
        default:
        case SPEX_MPZ:
            SPEX_CHECK(spex_demo_tripread_mpz(A_handle, file, option));
            break;
        case SPEX_FP64:
            SPEX_CHECK(spex_demo_tripread_double(A_handle, file, option));
            break;
    }

    return SPEX_OK;
}
