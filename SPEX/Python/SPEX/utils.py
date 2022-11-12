# SPEX: (c) 2022, Chris Lourenco, United States Naval Academy,
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

import numpy as np
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix, isspmatrix, isspmatrix_csc, linalg

def spex_matrix_from_file(fname):
    #fname is the name of the file that contains matrix A

    ##--------------------------------------------------------------------------
    ## Load file data and store it in mat
    ##--------------------------------------------------------------------------
    mat = np.loadtxt(fname, delimiter=' ')

    ##--------------------------------------------------------------------------
    ## Splice mat to get the componets of matrix A and populate a scipy sparse matrix
    ##--------------------------------------------------------------------------
    data = mat[1:,2]
    row = mat[1:,0].astype(int)
    col = mat[1:,1].astype(int)
    n = mat[0][0].astype(int)
    m = mat[0][1].astype(int)
    ## The matrix in the file is in triplet form
    triplet=coo_matrix((data, (row, col)),shape=(n, m))
    ## Change the triplet matrix into a csc matrix
    A = triplet.tocsc()

    return A

class Options:
    def __init__(self, out="double", ordering=None):
        self.ordering = ordering
        self.output = out

    def default_lu(self):
        self.ordering="colamd"

    def default_chol(self):
        self.ordering="amd"

    def order(self):

        if self.ordering=="none":
            order=0
        elif self.ordering=="colamd": ##colamd is the default ordering for Left LU
            order=1
        elif self.ordering=="amd": ##amd is the default ordering for Cholesky
            order=2
        else:
            print("Invalid order options")
            raise ValueError

        return order

    def charOut(self):

        if self.output=="double":
            charOut=False
        elif self.output=="string":
            charOut=True
        else:
            print("Invalid output type options")
            raise ValueError

        return charOut


