## demo a demo of SPEX_Chol
# SPEX_Chol is a package for solving sparse SPD linear systems of equations
# with a roundoff-free integer-preserving method.  The result is
# always exact, unless the matrix A is perfectly singular.

## SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
## Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
## Texas A&M University. All Rights Reserved. 
## SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#import SPEX_Chol
from SPEX import SPEX_Chol

#import scientific computing 
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix


def test():
    ##baby test with small dense matrix
    ##just to check the "connection" with c works
    row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
    data = np.array([4, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
    A=csc_matrix((data, (row, col)), shape=(3, 3))
    b=np.ones(3,dtype=np.float64)

    x=SPEX_Chol.Cholesky(A,b,{'SolutionType': 'string', 'Ordering':'amd'})
    for i in range(3):
        print(x[i])
    x=SPEX_Chol.Cholesky(A,b,{'SolutionType': 'double', 'Ordering':'amd'})
    for i in range(3):
        print(x[i])
        
    
    

##doesnt need to be a function, idk what is cleaner
def test2():
    ##--------------------------------------------------------------------------
    ## Read in A and populate b
    ##---------------------------------------------------------------------
    fname='../ExampleMats/1438.mat.txt'
    A=SPEX_Chol.spex_matrix_from_file(fname)
    b=np.ones(A.shape[0],dtype=np.float64)
    
    ##--------------------------------------------------------------------------
    ## solve
    ##--------------------------------------------------------------------------
    x=SPEX_Chol.Cholesky(A,b)
    print('Test from file successful')
    print(x)

test2()

test()

##doesnt need to be a function, idk what is cleaner
def test_void():
    ##--------------------------------------------------------------------------
    ## Read in A and populate b
    ##---------------------------------------------------------------------
    fname='../ExampleMats/1438.mat.txt'
    A=SPEX_Chol.spex_matrix_from_file(fname)
    b=np.ones(A.shape[0],dtype=np.float64)
    ##baby test with small dense matrix
    ##just to check the "connection" with c works
    row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
    data = np.array([4, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
    A=csc_matrix((data, (row, col)), shape=(3, 3))
    b=np.ones(3,dtype=np.float64)
    
    ##--------------------------------------------------------------------------
    ## solve
    ##--------------------------------------------------------------------------
    x=SPEX_Chol.SPEX_Chol_void(A,b,2)
    print('Test from file successful')
    for i in range(3):
        print(x[i])
        
test_void()