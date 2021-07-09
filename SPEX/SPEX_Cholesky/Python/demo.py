from SPEX_Choleskpy import SPEX_Chol_backslash

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix


##doesnt need to be a function, idk what is cleaner
def test2():
    ##--------------------------------------------------------------------------
    ## Read in A and populate b
    ##---------------------------------------------------------------------
    fname='../ExampleMats/1438.mat.txt'
    A=SPEX_Chol_backslash.spex_matrix_from_file(fname)
    b=np.ones(A.shape[0],dtype=np.float64)
    
    ##--------------------------------------------------------------------------
    ## solve
    ##--------------------------------------------------------------------------
    x=SPEX_Chol_backslash.SPEX_Chol_backslash(A,b)
    print('Test from file successful')

#test2()

SPEX_Chol_backslash.test()