//------------------------------------------------------------------------------
// Demos/SPEX_Backslash_demo: example of SPEX_Blackslash
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2021-2022, Chris Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* A demo of SPEX_Backslash in C
 */
 
# include "demos.h"

#define FREE_WORKSPACE                  \
{                                       \
    SPEX_matrix_free(&A,NULL);          \
    SPEX_matrix_free(&b,NULL);          \
    SPEX_matrix_free(&x,NULL);          \
    SPEX_FREE(option);                  \
    SPEX_finalize();                    \
}                                       \


int main( int argc, char *argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int64_t n = 0, ok;

    SPEX_matrix A = NULL;
    SPEX_matrix b = NULL;
    SPEX_matrix x = NULL;

    // Set default options
    SPEX_options option = NULL;
    OK (SPEX_create_default_options(&option));

    char *mat_name, *rhs_name;
    int64_t rat = 1;

    // Process the command line
    OK (SPEX_backslash_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    // Read in A
    FILE *mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }

    // Note, there are a few matrices in BasisLIB that dont fit in double
    // Need to use the other tripread for those.
    OK (SPEX_tripread_double(&A, mat_file, option));
    fclose(mat_file);
    n = A->n;

    FILE *rhs_file = fopen(rhs_name,"r");
    if( rhs_file == NULL )
    {
        printf("\nNo RHS file provided");
        printf("\nWill generate a RHS of all 1s");
        SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ, n, 1, n,
                             false, true, option);
        for (int64_t k = 0; k < A->n; k++)
            OK (SPEX_mpz_set_ui(b->x.mpz[k],1));

    }
    else
    {
        OK (SPEX_read_dense(&b, rhs_file, option));
        fclose(rhs_file);
    }


    //--------------------------------------------------------------------------
    // Solve Ax = b
    //--------------------------------------------------------------------------

    clock_t start = clock();

    option->print_level = 0;

    DEMO_OK( SPEX_backslash(&x, SPEX_FP64, A, b, option));

    // Chris:
    TRY ( SPEX_backslash(&x, SPEX_FP64, A, b, context));
    TRY ( SPEX_mpfr_fprintf(context, fp, "huh? %g\n", x))  ;

    // consistent:
    TRY ( SPEX_backslash(context, &x, SPEX_FP64, A, b));
    TRY ( SPEX_mpfr_fprintf(context, fp, ":'-| %g\n", x))  ;

    clock_t end = clock();

    double t_tot = (double) (end - start) / CLOCKS_PER_SEC;

    printf("\nSPEX Backslash Factor & Solve time: %lf\n", t_tot);

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
}

