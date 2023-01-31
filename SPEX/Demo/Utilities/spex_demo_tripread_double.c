//------------------------------------------------------------------------------
// Demo/SPEX_tripread_double: reads a double matrix stored in triplet format
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function reads a double matrix stored in triplet format
 * This format used can be seen in any of the example mat files.
 *
 * The first line of the file contains three integers: m, n, nnz,
 * where the matrix is m-by-n with nnz entries.
 *
 * This is followed by nnz lines, each containing a single triplet: i, j, aij,
 * which defines the row index (i), column index (j), and value (aij) of
 * the entry A(i,j).  The value aij is a floating-point number.
 */

#include "demos.h"

SPEX_info spex_demo_tripread_double
(
    SPEX_matrix *A_handle,      // Matrix to be populated
    FILE *file,                 // file to read from (must already be open)
    SPEX_options option
)
{

    SPEX_info info ;
    ASSERT(A_handle!=NULL);
    ASSERT(file!=NULL);
    
    (*A_handle) = NULL ;

    // Read in triplet form first
    int64_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    int s = fscanf(file, "%"PRId64" %"PRId64" %"PRId64"\n", &m, &n, &nz);
    if (feof(file) || s < 3)
    {
        printf ("premature end-of-file\n");
        return SPEX_INCORRECT_INPUT;
    }

    // First, we create our A matrix which is triplet double
    SPEX_matrix A = NULL;
    info = SPEX_matrix_allocate(&A, SPEX_TRIPLET, SPEX_FP64, m, n, nz,
        false, true, option);
    if (info != SPEX_OK)
    {
        printf ("unable to allocate matrix\n");
        return (info);
    }

//  s = fscanf (file, "%"PRId64" %"PRId64" %lf\n",
//      &(A->i[0]), &(A->j[0]), &(A->x.fp64[0]));
//  if (feof(file) || s < 2)
//  {
//      printf ("premature end-of-file\n");
//      SPEX_matrix_free(&A, option);
//      return SPEX_INCORRECT_INPUT;
//  }

//  // Matrices in this format are 1 based. We decrement
//  // the indices by 1 to use internally
//  A->i[0] -= 1;
//  A->j[0] -= 1;

    // Read in the values from file
    for (int64_t k = 0; k < nz; k++)
    {
        s = fscanf(file, "%"PRId64" %"PRId64" %lf\n",
            &(A->i[k]), &(A->j[k]), &(A->x.fp64[k]));
        if ((feof(file) && k != nz-1) || s < 3)
        {
            printf ("premature end-of-file\n");
            SPEX_matrix_free(&A, option);
            return SPEX_INCORRECT_INPUT;
        }
//      // Conversion from 1 based to 0 based
        A->i[k] -= 1;
        A->j[k] -= 1;
    }

    // the triplet matrix now has nz entries
    A->nz = nz;

// print_level from option struct:
//      0: nothing
//      1: just errors
//      2: errors and terse output
//      3: verbose

//  option->print_level = 3 ;
    info = SPEX_matrix_check (A, option);
    if (info != SPEX_OK)
    {
        printf ("invalid matrix\n");
        return (info);
    }

    // At this point, A is a double triplet matrix. We make a copy of it with C

    SPEX_matrix C = NULL;
    info = SPEX_matrix_copy(&C, SPEX_CSC, SPEX_MPZ, A, option);
    if (info != SPEX_OK)
    {
        printf ("unable to copy matrix\n");
        return (info);
    }

    // Success. Set A_handle = C and free A

    SPEX_matrix_free(&A, option);
    (*A_handle) = C;
    return (info);
}
