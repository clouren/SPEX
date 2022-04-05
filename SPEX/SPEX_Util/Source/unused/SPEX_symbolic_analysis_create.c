//------------------------------------------------------------------------------
// SPEX_Util/SPEX_symbolic_analysis_create
//------------------------------------------------------------------------------


// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function creates the SPEX_symbolic_analysis struct
 *
 * Input is the SPEX_symbolic_analysis structure, everything is NULL upon
 * function termination.
 */
//TODO coments


#include "spex_util_internal.h"

SPEX_info SPEX_symbolic_analysis_create
(
    SPEX_symbolic_analysis **S, // Structure to be created
    const SPEX_options *option
)
{
  if (!spex_initialized ( )) return (SPEX_PANIC) ;

  ((*S)->P_perm)=NULL;
  (*S)->Pinv_perm=NULL;
  (*S)->Q_perm=NULL;
  (*S)->Qinv_perm=NULL;

  (*S)->parent=NULL;
  ((*S)->cp)=NULL;
  ((*S)->c)=NULL;


  return SPEX_OK;
}