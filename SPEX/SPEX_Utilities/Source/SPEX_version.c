//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_version: report SPEX version information
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Returns the library version, date, and thread-safety status

// FIXME: add to user guide

#include "spex_util_internal.h"

SPEX_info SPEX_version
(
    // output
    int version [3],            // SPEX major, minor, and sub version
    char date [128],            // date of this version
    char thread_safety [128]    // a string describing the mechanism used to
                                // make SPEX thread-safe
)
{

    if (version != NULL)
    {
        version [0] = SPEX_VERSION_MAJOR ;
        version [1] = SPEX_VERSION_MINOR ;
        version [2] = SPEX_VERSION_SUB ;
    }

    if (date != NULL)
    {
        strncpy (date, SPEX_DATE, 127);
    }

    if (thread_safety != NULL)
    {
        #if defined ( SPEX_USE_PTHREADS )
            strncpy (thread_safety, "POSIX pthreads", 127);
        #elif defined ( SPEX_USE_WIN32_THREADS )
            strncpy (thread_safety, "Windows threads", 127);
        #elif defined ( _OPENMP )
            strncpy (thread_safety, "OpenMP", 127);
        #else
            strncpy (thread_safety, "unsafe", 127);
        #endif
    }

    return (SPEX_OK);
}

