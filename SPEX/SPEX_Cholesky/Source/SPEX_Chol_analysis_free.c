//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_analysis_free: Free SPEX_Chol_analysis data struct
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

/* Purpose: Free the SPEX_Chol_analysis structure*/
void SPEX_Chol_analysis_free
(
    SPEX_Chol_analysis** S
)
{
  if ((S != NULL) && (*S != NULL))
  {
    SPEX_FREE ((*S)->p) ;
    SPEX_FREE((*S)->parent);
    SPEX_FREE((*S)->cp);
    SPEX_FREE((*S)->pinv);
    SPEX_FREE (*S) ;
  }
}
