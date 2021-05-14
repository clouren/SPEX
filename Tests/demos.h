//------------------------------------------------------------------------------
// SPEX_Chol/Demo/demos.h: #include file the demo programs
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_Util.h"
#include "SPEX_Left_LU.h"
#include "SPEX_Chol.h"
#include "SPEX_WAMF.h"
#include <time.h>
#include <stdint.h>
#include <inttypes.h>

#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))
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
void SPEX_Chol_determine_error
(
    SPEX_info ok
);

SPEX_info SPEX_tripread
(
    SPEX_matrix **A_handle,      // Matrix to be constructed
    FILE* file,                  // file to read from (must already be open)
    SPEX_options* option         // Command options
);