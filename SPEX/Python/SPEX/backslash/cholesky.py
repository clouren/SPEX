#-------------------------------------------------------------------------------
# SPEX/Python/backslash/cholesky.py: solve Ax=b using Cholesky factorization
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

def cholesky( A, b, options=Options('double', 'amd')):
    ## A is a scipy.sparse(data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    ## options is a dictionary that specifies what tipe the solution should be, this by default is double

    ##--------------------------------------------------------------------------
    ## Verify inputs
    ##--------------------------------------------------------------------------
    if not isspmatrix(A):
        raise TypeError("Input matrix must be scipy.sparse")
    ## If the sparse input matrix is not in csc form, convert it into csc form
    if not isspmatrix_csc(A):
        A.tocsc()
    ## Check symmetry
    tol=1e-8
    if scipy.sparse.linalg.norm(A-A.T, scipy.Inf) > tol:
        raise TypeError("Input matrix is not symmetric")
    # Check input shape
    if A.shape[1]!=b.shape[0]:
        raise TypeError("Matrix input shapes must match")

    if options.ordering==None:
        options.default_chol()

    ##--------------------------------------------------------------------------
    ## Call SPEX
    ##--------------------------------------------------------------------------
    x=spex_connect(A,b,options.order(),options.charOut(),3) #3 calls the Cholesky

    return x
