#include "spex_left_lu_internal.h"

SPEX_info SPEX_LU_analyze
(
    SPEX_symbolic_analysis** S_handle, // symbolic analysis including
                                 // column perm. and nnz of L and U
    const SPEX_matrix *A,        // Input matrix
    const SPEX_options *option   // Control parameters, if NULL, use default
)
{
    return spex_lu_analyze(S_handle, A, option);
}
