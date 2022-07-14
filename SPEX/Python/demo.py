## demo.py a demo of SPEX_Chol

# SPEX: (c) 2022, Chris Lourenco, United States Naval Academy, 
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved. 
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

# SPEX is a package for solving sparse linear systems of equations
# with a roundoff-free integer-preserving method.  The result is
# always exact, unless the matrix A is perfectly singular.

# Import SPEX
from SPEX import utils, backslash
from SPEX.utils import Options

# Import scientific computing
import numpy as np
from numpy.random import default_rng
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import random
from scipy import stats


##--------------------------------------------------------------------------
## Cholesky
##--------------------------------------------------------------------------

# Read in A from file and populate b
fname='../SPEX_Cholesky/ExampleMats/1438.mat.txt'
A=utils.spex_matrix_from_file(fname)
b=np.ones(A.shape[0],dtype=np.float64)
# Solve
x=backslash.cholesky(A,b)
print(x)


##--------------------------------------------------------------------------
## Left LU
##--------------------------------------------------------------------------
# Create A and B
row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
data = np.array([4, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
A=csc_matrix((data, (row, col)), shape=(3, 3))
b=np.ones(3,dtype=np.float64)

# Solve
options=Options("string")
x=backslash.left_lu(A,b,options)
print(x)

##--------------------------------------------------------------------------
## Backslash
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
x=backslash.general(A,b)
print(x)