//------------------------------------------------------------------------------
// SPEX_Backslash/Demo/demos.c: support functions for the demo programs
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2021, Chris Lourenco, Lorena Mejia Domenzain, Erick
// Moreno-Centeno, Timothy A. Davis, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_backslash_process_command_line: process command line for demo programs.
// SPEX_tripread_double: read a double matrix from a file in triplet format.
// SPEX_read_dense: read a dense matrix from a file.

// FIXME: duplicate functions.  See lu_demos.h
// FIXME: one function per file

#include "back_demos.h"

// ignore warnings about unused parameters in this file
#pragma GCC diagnostic ignored "-Wunused-parameter"

//------------------------------------------------------------------------------
// SPEX_backslash_process_command_line
//------------------------------------------------------------------------------

/* Purpose: This processes the command line for user specified options */
SPEX_info SPEX_backslash_process_command_line //processes the command line
(
    int64_t argc,           // number of command line arguments
    char *argv[],           // set of command line arguments
    SPEX_options option,    // struct containing the command options
    char **mat_name,        // Name of the matrix to be read in
    char **rhs_name,        // Name of the RHS vector to be read in
    int64_t *rat            // data type of output solution.
                            // 1: mpz, 2: double, 3: mpfr
)
{
    //SPEX_info ok;
    for (int64_t  i = 1; i < argc; i++)
    {
        char *arg = (char*) argv[i];
        if ( strcmp(arg,"help") == 0)
        {
            //SPEX_show_usage();
            return SPEX_INCORRECT_INPUT;
        }
        // Type of pivoting if LU is used
        else if ( strcmp(arg,"p") == 0 || strcmp(arg,"piv") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be a pivot argument between"
                    " 0-5 following p\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->pivot = atoi(argv[i]);
            if (option->pivot < 0 || option->pivot > 5)
            {
                printf("\n****ERROR! Invalid pivot selection!"
                    "\nDefaulting to smallest pivot\n\n");
                option->pivot = SPEX_SMALLEST;
            }
        }
        // Type of column ordering used for LU or row/col permutation for Chol
        else if ( strcmp(arg, "q") == 0 || strcmp(arg,"col") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be an argument between 0-2"
                    "following q\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->order = atoi(argv[i]);
            if (option->order < 0 || option->order > 2)
            {
                printf("\n****ERROR! Invalid column orderint64_t*g"
                    "\nDefaulting to COLAMD\n\n");
                option->order = SPEX_COLAMD;
            }
        }
        // Tolerance for pivoting in LU
        else if ( strcmp(arg,"t") == 0 || strcmp(arg, "tol") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be a non negative tolerance"
                    " value following t\n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atof(argv[i]))
            {
                printf("\n****ERROR! There must be a non negative tolerance"
                    " value following t\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->tol = atof(argv[i]);
            if (option->tol < 0)
            {
                printf("\n****ERROR! Invalid Tolerance, tolerance must be"
                    " non-negative\n");
                return SPEX_INCORRECT_INPUT;
            }
        }
        else if ( strcmp(arg,"out2") == 0 || strcmp(arg, "o2") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print nothing) 1 (print err) or 2 (terse) \n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (printt nothing) 1 (print err) or 2 (terse) \n");
                return SPEX_INCORRECT_INPUT;
            }
            option->print_level = atoi(argv[i]);
        }
        else if ( strcmp(arg, "out") == 0 || strcmp(arg, "o") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o or out must be followed by"
                    " 1 (rational) 2 (double) or 3 (variable precision) \n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o or out must be followed by"
                    " 1 (rational) 2 (double) or 3 (variable precision) \n");
                return SPEX_INCORRECT_INPUT;
            }
            *rat = atoi(argv[i]);
            if (*rat < 1 || *rat > 3)
            {
                printf("\n\n****ERROR! Invalid output type!\n"
                   "Defaulting to rational");
                *rat = 1;
            }
            if (*rat == 3)
            {
                if (!argv[++i])
                {
                    printf("\n****ERROR! Precision must be specified\n");
                    return SPEX_INCORRECT_INPUT;
                }
                else if (!atoi(argv[i]))
                {
                    printf("\n****ERROR! Precision must be specified\n");
                    return SPEX_INCORRECT_INPUT;
                }
                option->prec = atoi(argv[i]);
                if (option->prec < 2)
                {
                    printf("\n\n****ERROR! Invalid precision. prec >= 2\n");
                    return SPEX_INCORRECT_INPUT;
                }
            }
        }
        else if ( strcmp(arg, "f") == 0 || strcmp(arg, "file") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! Matrix name must be entered\n");
                return SPEX_INCORRECT_INPUT;
            }
            *mat_name = argv[i];
            if (!argv[++i])
            {
                printf("\n****ERROR! Right hand side vector name must"
                    " be entered\n");
                return SPEX_INCORRECT_INPUT;
            }
            *rhs_name = argv[i];
        }
        else
        {
            printf("\n\n**ERROR! Unknown command line parameter: %s"
                    "\nIgnoring this parameter\n",arg);
            return SPEX_INCORRECT_INPUT;
        }
    }
    return SPEX_OK;
}

//------------------------------------------------------------------------------
// SPEX_tripread
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

SPEX_info SPEX_tripread
(
    SPEX_matrix *A_handle,      // Matrix to be constructed
    FILE *file,                 // file to read from (must already be open)
    SPEX_options option         // Command options
)
{

    SPEX_info info ;
    if (A_handle == NULL || file == NULL)
    {
        printf ("invalid input\n") ;
        return SPEX_INCORRECT_INPUT;
    }

    (*A_handle) = NULL ;

    int64_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    int s = fscanf(file, "%"PRId64" %"PRId64" %"PRId64"\n", &m, &n, &nz);
    if (feof(file) || s < 3)
    {
        printf ("premature end-of-file\n") ;
        return SPEX_INCORRECT_INPUT;
    }

    // Allocate memory for A
    // A is a triplet mpz_t matrix
    SPEX_matrix A = NULL;
    info = SPEX_matrix_allocate(&A, SPEX_TRIPLET, SPEX_MPZ, m, n, nz,
        false, true, option);
    if (info != SPEX_OK)
    {
        return (info) ;
    }

    // Read in first values of A
    info = SPEX_gmp_fscanf(file, "%"PRId64" %"PRId64" %Zd\n",
        &A->i[0], &A->j[0], &A->x.mpz[0]);
    if (feof (file) || info != SPEX_OK)
    {
        printf ("premature end-of-file\n") ;
        SPEX_matrix_free(&A, option);
        return SPEX_INCORRECT_INPUT;
    }

    // Matrices in this format are 1 based, so we decrement by 1 to get
    // 0 based for internal functions
    A->i[0] -= 1;
    A->j[0] -= 1;

    // Read in the values from file
    for (int64_t p = 1; p < nz; p++)
    {
        info = SPEX_gmp_fscanf(file, "%"PRId64" %"PRId64" %Zd\n",
            &A->i[p], &A->j[p], &A->x.mpz[p]);
        if ((feof(file) && p != nz-1) || info != SPEX_OK)
        {
            printf ("premature end-of-file\n") ;
            SPEX_matrix_free(&A, option);
            return SPEX_INCORRECT_INPUT;
        }
        // Conversion from 1 based to 0 based if necessary
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
    return (info) ;
}

//------------------------------------------------------------------------------
// SPEX_tripread_double
//------------------------------------------------------------------------------

/* Purpose: This function reads in a matrix stored in a triplet format
 * with double entries. The format used can be seen in any of the
 * example mat files.
 *
 * This is only used for Demo purposes
 */

SPEX_info SPEX_tripread_double
(
    SPEX_matrix *A_handle,      // Matrix to be populated
    FILE *file,                 // file to read from (must already be open)
    SPEX_options option         // Command options
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

    if (feof(file) || s <= 0)
    {
        printf ("premature end-of-file\n") ;
        SPEX_matrix_free(&A, option);
        return SPEX_INCORRECT_INPUT;
    }

    // Matrices in this format are either 1 or 0 based.
    // We check the decrement for the indices
    int64_t decrement;
    if (SPEX_MIN(A->i[0], A->j[0]) == 0)
        decrement = 0;
    else
        decrement = 1;
    A->i[0] -= decrement;
    A->j[0] -= decrement;

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
        // Conversion from 1 based to 0 based if needed
        A->i[k] -= decrement;
        A->j[k] -= decrement;
    }

    // the triplet matrix now has nz entries
    A->nz = nz;

    // At this point, A is a double triplet matrix. We make a copy of it with C
    // C is a CSC matrix with mpz entries

    SPEX_matrix C = NULL;
    SPEX_matrix_copy(&C, SPEX_CSC, SPEX_MPZ, A, option);

    // Success. Set A_handle = C and free A

    SPEX_matrix_free(&A, option);
    (*A_handle) = C;
    return (info) ;
}

//------------------------------------------------------------------------------
// SPEX_read_dense
//------------------------------------------------------------------------------

/* Purpose: Read a dense matrix for RHS vectors.
 * the values in the file must be integers
 */

SPEX_info SPEX_read_dense
(
    SPEX_matrix *b_handle,  // Matrix to be constructed
    FILE *file,             // file to read from (must already be open)
    SPEX_options option
)
{

    if (file == NULL)
    {
        printf ("invalid inputs\n") ;
        return SPEX_INCORRECT_INPUT;
    }
    int64_t nrows, ncols;
    SPEX_info info ;

    // First, we obtain the dimension of the matrix
    int s = fscanf(file, "%"PRId64" %"PRId64, &nrows, &ncols) ;
    if (feof(file) || s < 2)
    {
        printf ("premature end-of-file\n") ;
        return SPEX_INCORRECT_INPUT;
    }

    // Now, we create our dense mpz_t matrix
    SPEX_matrix A = NULL;
    info = SPEX_matrix_allocate(&A, SPEX_DENSE, SPEX_MPZ, nrows, ncols,
        nrows*ncols, false, true, option);
    if (info != SPEX_OK)
    {
        return (info) ;
    }

    // We now populate the matrix b.
    for (int64_t i = 0; i < nrows; i++)
    {
        for (int64_t j = 0; j < ncols; j++)
        {
            info = SPEX_gmp_fscanf(file, "%Zd", &(SPEX_2D(A, i, j, mpz)));
            if (info != SPEX_OK)
            {
                printf("\n\nhere at i = %"PRId64" and j = %"PRId64"", i, j);
                return SPEX_INCORRECT_INPUT;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Success, set b_handle = A
    //--------------------------------------------------------------------------

    (*b_handle) = A;
    return (info) ;
}
