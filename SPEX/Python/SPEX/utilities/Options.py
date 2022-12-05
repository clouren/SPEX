# SPEX: (c) 2022, Chris Lourenco, United States Naval Academy,
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

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


