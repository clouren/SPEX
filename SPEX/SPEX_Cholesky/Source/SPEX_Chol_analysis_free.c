//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_analysis_free: Free SPEX_Chol_analysis data struct
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"
/* Purpose: Free the SPEX_Chol_analysis Structure*/
void SPEX_Chol_analysis_free
(
    SPEX_Chol_analysis* S
)
{
    SPEX_FREE(S->parent);
    SPEX_FREE(S->cp);
    SPEX_FREE(S);
}
