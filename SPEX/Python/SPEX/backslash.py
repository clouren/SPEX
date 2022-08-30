# SPEX: (c) 2022, Chris Lourenco, United States Naval Academy, 
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved. 
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

from SPEX.utils import Options
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

def left_lu( A, b, options=Options('double', 'colamd')): 
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
        
    if options.ordering==None:
        options.default_lu()
        
    ##--------------------------------------------------------------------------
    ## Call SPEX
    ##--------------------------------------------------------------------------
    x=spex_connect(A,b,options.order(),options.charOut(),2) #3 calls lu

    return x

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
