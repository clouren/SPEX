//------------------------------------------------------------------------------
// SPEX_tripread_double
//------------------------------------------------------------------------------

/* Purpose: This function reads in a double matrix stored in a triplet format
 * This format used can be seen in any of the example mat files.
 *
 * The first line of the file contains three integers: m, n, nnz,
 * where the matrix is m-by-n with nnz entries.
 *
 * This is followed by nnz lines, each containing a single triplet: i, j, aij,
 * which defines the row index (i), column index (j), and value (aij) of
 * the entry A(i,j).  The value aij is a floating-point number.
 */

#include "SPEX.h"

SPEX_info SPEX_tripread_double
(
    SPEX_matrix *A_handle,      // Matrix to be populated
    FILE *file,                 // file to read from (must already be open)
    SPEX_options option
)
{

    SPEX_info info ;
    if (A_handle == NULL || file == NULL)
    {
        printf ("invalid input\n") ;
        return SPEX_INCORRECT_INPUT;
    }
    (*A_handle) = NULL ;

    // Read in triplet form first
    int64_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    int s = fscanf(file, "%"PRId64" %"PRId64" %"PRId64"\n", &m, &n, &nz);
    if (feof(file) || s < 3)
    {
        printf ("premature end-of-file\n") ;
        return SPEX_INCORRECT_INPUT;
    }

    // First, we create our A matrix which is triplet double
    SPEX_matrix A = NULL;
    info = SPEX_matrix_allocate(&A, SPEX_TRIPLET, SPEX_FP64, m, n, nz,
        false, true, option);
    if (info != SPEX_OK)
    {
        return (info) ;
    }

    s = fscanf (file, "%"PRId64" %"PRId64" %lf\n",
        &(A->i[0]), &(A->j[0]), &(A->x.fp64[0])) ;
    if (feof(file) || s < 2)
    {
        printf ("premature end-of-file\n") ;
        SPEX_matrix_free(&A, option);
        return SPEX_INCORRECT_INPUT;
    }

    // Matrices in this format are 1 based. We decrement
    // the indices by 1 to use internally
    A->i[0] -= 1;
    A->j[0] -= 1;

    // Read in the values from file
    for (int64_t k = 1; k < nz; k++)
    {
        s = fscanf(file, "%"PRId64" %"PRId64" %lf\n",
            &(A->i[k]), &(A->j[k]), &(A->x.fp64[k]));
        if ((feof(file) && k != nz-1) || s < 3)
        {
            printf ("premature end-of-file\n") ;
            SPEX_matrix_free(&A, option);
            return SPEX_INCORRECT_INPUT;
        }
        // Conversion from 1 based to 0 based
        A->i[k] -= 1;
        A->j[k] -= 1;
    }

    // the triplet matrix now has nz entries
    A->nz = nz;

    // At this point, A is a double triplet matrix. We make a copy of it with C

    SPEX_matrix C = NULL;
    SPEX_matrix_copy(&C, SPEX_CSC, SPEX_MPZ, A, option);

    // Success. Set A_handle = C and free A

    SPEX_matrix_free(&A, option);
    (*A_handle) = C;
    return (info) ;
}
