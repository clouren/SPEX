//------------------------------------------------------------------------------
// SPEX_Util/SPEX_factorization_create
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

//TODO coments fix

#include "spex_util_internal.h"

SPEX_info SPEX_factorization_create
(
    SPEX_factorization **F, // Structure to be created
    const SPEX_options *option
)
{

  if (!spex_initialized ( )) return (SPEX_PANIC) ;

  (*F)->L=NULL;
  (*F)->U=NULL;
  (*F)->rhos=NULL;

  (*F)->Q=NULL;
  (*F)->R=NULL;

  (*F)->P_perm= NULL;
  (*F)->Pinv_perm= NULL;
  (*F)->Q_perm= NULL;
  (*F)->Qinv_perm= NULL;


  return SPEX_OK;
}