#include "spex_left_lu_internal.h"

SPEX_info SPEX_Left_LU_factorize
(
    // output:
    SPEX_factorization **F_handle, // LU factorization
    // input:
    const SPEX_matrix *A,      // matrix to be factored
    const SPEX_symbolic_analysis *S, // symbolic analysis
    const SPEX_options* option // command options
)
{
    return spex_left_lu_factorize(F_handle, A, S, option);
}
