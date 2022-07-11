//------------------------------------------------------------------------------
// SPEX_Chol/Demo/demos.c: support functions for the demo programs
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_Chol_process_command_line: process command line for demo programs.
// SPEX_tripread_double: read a double matrix from a file in triplet format.
// SPEX_read_dense: read a dense matrix from a file.
// SPEX_Chol_determine_error: error codes for exact factorization

#include "demos.h"
#include "spex_gmp.h"
#include "spex_util_internal.h"


// ignore warnings about unused parameters in this file
#pragma GCC diagnostic ignored "-Wunused-parameter"

//------------------------------------------------------------------------------
// SPEX_Chol_process_command_line
//------------------------------------------------------------------------------

/* Purpose: This processes the command line for user specified options */
SPEX_info SPEX_Chol_process_command_line //processes the command line
(
    int64_t argc,           // number of command line arguments
    char* argv[],           // set of command line arguments
    SPEX_options* option,   // struct containing the command options
    char** mat_name,        // Name of the matrix to be read in
    char** rhs_name,        // Name of the RHS vector to be read in
    int64_t *rat            // data type of output solution.
                            // 1: mpz, 2: double, 3: mpfr
)
{
    //SPEX_info ok;
    for (int64_t  i = 1; i < argc; i++)
    {
        char* arg = (char*) argv[i];
        if ( strcmp(arg,"help") == 0)
        {
            //SPEX_show_usage();
            return SPEX_INCORRECT_INPUT;
        }
        else if ( strcmp(arg, "q") == 0 || strcmp(arg,"col") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be an argument between 0-2"
                    "followint64_t*g q\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->order = atoi(argv[i]);
            if (option->order < 0 || option->order > 2)
            {
                printf("\n****ERROR! Invalid column orderint64_t*g"
                    "\nDefaultint64_t*g to COLAMD\n\n");
                option->order = SPEX_COLAMD;
            }
        }
        else if ( strcmp(arg,"out2") == 0 || strcmp(arg, "o2") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print64_t*64_tnothint64_t*g) 1 (print64_t*64_terr) or 2 (terse) \n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print64_t*64_tnothint64_t*g) 1 (print64_t*64_terr) or 2 (terse) \n");
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
                   "Defaultint64_t*g to rational");
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
            printf("\n\n**ERROR! Unknown command lint64_t*e parameter: %s"
                    "\nIgnorint64_t*g this parameter\n",arg);
            return SPEX_INCORRECT_INPUT;
        }
    }
    return SPEX_OK;
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
    SPEX_matrix **A_handle,     // Matrix to be populated
    FILE* file,                 // file to read from (must already be open)
    SPEX_options* option        // Command options
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
    SPEX_matrix *A = NULL;
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
    // C is a CSC matrix with mpz entries
    
    SPEX_matrix* C = NULL;
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
    SPEX_matrix **b_handle, // Matrix to be constructed
    FILE* file,             // file to read from (must already be open)
    SPEX_options* option
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
    SPEX_matrix* A = NULL;
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
            info = SPEX_gmp_fscanf(file, "%Zd", &(A->x.mpz[i+j*nrows]));
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

/* Purpose: Determine why a SPEX_Chol function failed
 */
void SPEX_Chol_determine_error
(
    SPEX_info ok
)
{
    if (ok == SPEX_OUT_OF_MEMORY)
    {
        printf("\nOut of memory\n");
    }
    else if (ok == SPEX_SINGULAR)
    {
        printf("\nInput matrix is singular OR no diagonal pivot. Please ensure input is SPD\n");
    }
    else if (ok == SPEX_UNSYMMETRIC)
    {
        printf("\nInput matrix is unsymmetric, please try left LU\n");
    }
    else if (ok == SPEX_INCORRECT_INPUT)
    {
        printf("\nIncorrect input for a SPEX_Chol Function\n");
    }
}



SPEX_info SPEX_check_solution
(
    const SPEX_matrix *A,         // Input matrix
    const SPEX_matrix *x,         // Solution vectors
    const SPEX_matrix *b,         // Right hand side vectors
    const SPEX_options* option    // Command options
)
{
    if (!SPEX_initialize ( )) return (SPEX_PANIC) ;

    //--------------------------------------------------------------------------
    // check inputs. Input are also checked by the two callers
    //--------------------------------------------------------------------------

    SPEX_info info ;
    /*SPEX_REQUIRE (A, SPEX_CSC,   SPEX_MPZ) ;
    SPEX_REQUIRE (x, SPEX_DENSE, SPEX_MPQ) ;
    SPEX_REQUIRE (b, SPEX_DENSE, SPEX_MPZ) ;*/

    //--------------------------------------------------------------------------
    // Declare vars
    //--------------------------------------------------------------------------

    int64_t p, j, i, nz;
    SPEX_matrix *b2 = NULL;   // b2 stores the solution of A*x
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);
    mpq_t scale; SPEX_MPQ_SET_NULL(scale);

    SPEX_CHECK (SPEX_mpq_init(temp));
    SPEX_CHECK(SPEX_mpq_init(scale));
    SPEX_CHECK (SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPQ, b->m, b->n,
        b->nzmax, false, true, option));


    //--------------------------------------------------------------------------
    // perform SPEX_mpq_addmul in loops to compute b2 = A'*x, where A' is the
    // scaled matrix with all entries in integer
    //--------------------------------------------------------------------------

    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            for (p = A->p[i]; p < A->p[i + 1]; p++)
            {
                // temp = A[p][i] (note this must be done seperately since A is
                // mpz and temp is mpq)
                SPEX_CHECK(SPEX_mpq_set_z(temp, A->x.mpz[p]));

                // temp = temp*x[i]
                SPEX_CHECK(SPEX_mpq_mul(temp, temp,
                                        SPEX_2D(x, i, j, mpq)));

                // b2[p] = b2[p]-temp
                SPEX_CHECK(SPEX_mpq_add(SPEX_2D(b2, A->i[p], j, mpq),
                                        SPEX_2D(b2, A->i[p], j, mpq),temp));
            }
        }
    }

    //--------------------------------------------------------------------------
    // Apply scales of A and b to b2 before comparing the b2 with scaled b'
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_mpq_div(scale, b->scale, A->scale));

    // Apply scaling factor, but ONLY if it is not 1
    int r;
    SPEX_CHECK(SPEX_mpq_cmp_ui(&r, scale, 1, 1));
    if (r != 0)
    {
        nz = x->m * x->n;
        for (i = 0; i < nz; i++)
        {
            SPEX_CHECK(SPEX_mpq_mul(b2->x.mpq[i], b2->x.mpq[i], scale));
        }
    }

    //--------------------------------------------------------------------------
    // check if b==b2
    //--------------------------------------------------------------------------

    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            // temp = b[i] (correct b)
            SPEX_CHECK(SPEX_mpq_set_z(temp, SPEX_2D(b, i, j, mpz)));

            // set check false if b!=b2
            SPEX_CHECK(SPEX_mpq_equal(&r, temp, SPEX_2D(b2, i, j, mpq)));
            if (r == 0)
            {
                //printf("ERROR\n");
                info = SPEX_PANIC;
                j = b->n;
                break;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Print info
    //--------------------------------------------------------------------------

    int pr = SPEX_OPTION_PRINT_LEVEL (option) ;
    if (info == SPEX_OK)
    {
        SPEX_PR1 ("Solution is verified to be exact.\n") ;
    }
    else if (info == SPEX_PANIC)
    {
        // This can never happen.
        SPEX_PR1 ("ERROR! Solution is wrong. This is a bug; please "
                  "contact the authors of SPEX.\n") ;
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SPEX_MPQ_CLEAR(temp);                   \
    SPEX_MPQ_CLEAR(scale);                  \
    SPEX_matrix_free(&b2, NULL);
    return info;
}

