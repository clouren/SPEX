//------------------------------------------------------------------------------
// SPEX_Chol/Demo/demos.h: #include file the demo programs
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy,
// Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "SPEX.h"

#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))
/* Purpose: This processes the command line for user specified options */
SPEX_info SPEX_cholesky_process_command_line //processes the command line
(
    int64_t argc,           // number of command line arguments
    char* argv[],           // set of command line arguments
    SPEX_options* option,   // struct containing the command options
    char** mat_name,        // Name of the matrix to be read in
    char** rhs_name,        // Name of the RHS vector to be read in
    int64_t *rat            // data type of output solution.
                            // 1: mpz, 2: double, 3: mpfr
);

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
);

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
);

//------------------------------------------------------------------------------
// SPEX_read_dense
//------------------------------------------------------------------------------

/* Purpose: Print error codes for Chol factorization
 */
void SPEX_cholesky_determine_error
(
    SPEX_info ok
);


SPEX_info SPEX_check_solution
(
    const SPEX_matrix *A,         // Input matrix
    const SPEX_matrix *x,         // Solution vectors
    const SPEX_matrix *b,         // Right hand side vectors
    const SPEX_options* option    // Command options
);
