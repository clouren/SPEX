#-------------------------------------------------------------------------------
# SPEX/Python/spex_python_demo.py: demo of 4 backslash functions with different input
#                           matrices
#-------------------------------------------------------------------------------

# SPEX: (c) 2022-2025, Christopher Lourenco, Jinhao Chen,
# Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#------------------------------------------------------------------------------

# SPEX is a package for solving sparse linear systems of equations
# with a roundoff-free integer-preserving method.  The result is
# always exact, unless the matrix A is perfectly singular.

# Import SPEX
import SPEXpy as SPEX
from SPEXpy import Options

# Import scientific computing
import numpy as np
from numpy.random import default_rng
from scipy.sparse import csc_matrix
from scipy.sparse import random
from scipy import stats


##--------------------------------------------------------------------------
## Cholesky
##--------------------------------------------------------------------------

# Create A and B
row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
data = np.array([4, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
A=csc_matrix((data, (row, col)), shape=(3, 3))
b=np.ones(3,dtype=np.float64)

# Solve
print("\nTesting Python Cholesky\n")
x=SPEX.cholesky_backslash(A,b)
print(x)
print("\nCholesky success")
##-------------------------------------------------------------------------
## LDL
##------------------------------------------------------------------------
# Create A and B
row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
data = np.array([-3, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
A=csc_matrix((data, (row, col)), shape=(3, 3))
b=np.ones(3,dtype=np.float64)

# Solve
print("\nTesting Python LDL\n")
x=SPEX.ldl_backslash(A,b)
print(x)
print("\nLDL success")

##--------------------------------------------------------------------------
## Left LU
##--------------------------------------------------------------------------

# Generate a random sparse matrix A and populate b
n=10
rng = default_rng()
rvs = stats.poisson(25, loc=10).rvs
S = random(n, n, density=0.7, random_state=rng, data_rvs=rvs)
S2=S+np.eye(n)
A=csc_matrix(S2)
b=np.ones(n,dtype=np.float64)

# Solve
print("\nTesting Python LU\n")
options=Options("string")
x=SPEX.lu_backslash(A,b,options)
print(x)
print("\nLU success")

##--------------------------------------------------------------------------
## Backslash
##--------------------------------------------------------------------------

# Use the previous matrices

# Solve
print("\nTesting Python Backslash\n")
x=SPEX.backslash(A,b)
print(x)
print("\nBackslash success\n")
#SPEX.backslash always returns the output as float64
