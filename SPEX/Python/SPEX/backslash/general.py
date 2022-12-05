#-------------------------------------------------------------------------------
# SPEX/Python/backslash/general.py: solve Ax=b
#-------------------------------------------------------------------------------

# SPEX: (c) 2022, Chris Lourenco, Jinhao Chen,
# Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#------------------------------------------------------------------------------

from SPEX.utilities import Options
from SPEX._connect import spex_connect

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix, isspmatrix, isspmatrix_csc, linalg

def general( A, b, options=Options('double')):
    ## A is a scipy.sparse(data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    ## options is a dictionary that specifies what tipe the solution should be, this by default is double

    ##--------------------------------------------------------------------------
    ## Verify inputs
    ##--------------------------------------------------------------------------
    if not isspmatrix(A):
        print("Input matrix must be scipy.sparse")
        raise TypeError
    ## If the sparse input matrix is not in csc form, convert it into csc form
    if not isspmatrix_csc(A):
        A.tocsc()
    # Check input shape
    if A.shape[1]!=b.shape[0]:
        print("Matrix input shapes must match")
        raise TypeError

    ##--------------------------------------------------------------------------
    ## Call SPEX
    ##--------------------------------------------------------------------------
    x=spex_connect(A,b,0,options.charOut(),1) #3 calls the general Backslash

    return x

