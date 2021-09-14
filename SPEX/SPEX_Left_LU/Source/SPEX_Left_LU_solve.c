#include "spex_left_lu_internal.h"

SPEX_info SPEX_Left_LU_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_matrix **x_handle,  // rational solution to the system
    // input:
    const SPEX_matrix *b,   // right hand side vector
    const SPEX_factorization* F, // LU factorization
    const SPEX_options* option // Command options
)
{
    return spex_left_lu_solve(x_handle, b, F, option);
}
