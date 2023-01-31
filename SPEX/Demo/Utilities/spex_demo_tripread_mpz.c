//------------------------------------------------------------------------------
// Demo/SPEX_tripread
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function reads in a matrix stored in a triplet format
 * This format used can be seen in any of the example mat files.
 *
 * The first line of the file contains three integers: m, n, nnz,
 * where the matrix is m-by-n with nnz entries.
 *
 * This is followed by nnz lines, each containing a single triplet: i, j, aij,
 * which defines the row index (i), column index (j), and value (aij) of
 * the entry A(i,j).  The value aij is an integer.
 */

#include "demos.h"


SPEX_info spex_demo_tripread_mpz
(
    SPEX_matrix *A_handle,      // Matrix to be constructed
    FILE *file,                 // file to read from (must already be open)
    SPEX_options option         // Command options
)
{

    SPEX_info info ;
    ASSERT(A_handle!=NULL);
    ASSERT(file!=NULL);

    (*A_handle) = NULL ;

    int64_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    int s = fscanf(file, "%"PRId64" %"PRId64" %"PRId64"\n", &m, &n, &nz);
    if (feof(file) || s < 3)
    {
        printf ("premature end-of-file\n");
        return SPEX_INCORRECT_INPUT;
    }

    // Allocate memory for A
    // A is a triplet mpz_t matrix
    SPEX_matrix A = NULL;
    info = SPEX_matrix_allocate(&A, SPEX_TRIPLET, SPEX_MPZ, m, n, nz,
        false, true, option);
    if (info != SPEX_OK)
    {
        return (info);
    }

//  // Read in first values of A
//  info = SPEX_gmp_fscanf(file, "%"PRId64" %"PRId64" %Zd\n",
//      &A->i[0], &A->j[0], &A->x.mpz[0]);
//  if (feof (file) || info != SPEX_OK)
//  {
//      printf ("premature end-of-file\n");
//      SPEX_matrix_free(&A, option);
//      return SPEX_INCORRECT_INPUT;
//  }

//  // Matrices in this format are 1 based, so we decrement by 1 to get
//  // 0 based for internal functions
//  A->i[0] -= 1;
//  A->j[0] -= 1;

    // Read in the values from file
    for (int64_t p = 0; p < nz; p++)
    {
        info = SPEX_gmp_fscanf(file, "%"PRId64" %"PRId64" %Zd\n",
            &A->i[p], &A->j[p], &A->x.mpz[p]);
        if ((feof(file) && p != nz-1) || info != SPEX_OK)
        {
            printf ("premature end-of-file\n");
            SPEX_matrix_free(&A, option);
            return SPEX_INCORRECT_INPUT;
        }
//      // Conversion from 1 based to 0 based if necessary
        A->i[p] -= 1;
        A->j[p] -= 1;
    }

    // the triplet matrix now has nz entries
    A->nz = nz;

    // A now contains our input matrix in triplet format. We now
    // do a matrix copy to get it into CSC form
    // C is a copy of A which is CSC and mpz_t
    SPEX_matrix C = NULL;
    SPEX_matrix_copy(&C, SPEX_CSC, SPEX_MPZ, A, option);

    // Free A, set A_handle
    SPEX_matrix_free(&A, option);
    (*A_handle) = C;
    return (info);
}
