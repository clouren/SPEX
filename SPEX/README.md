SPEX is a software package for SParse EXact algebra

Files and folders in this distribution:

    README.md       this file

    build           Contains the SPEX C library as well
                    as the .so files

    Doc             User guide for the SPEX software package

    Demo            demo programs for SPEX

    ExampleMats     test matrices for demo programs

    Include         the SPEX.h include file for user applications

    MATLAB          MATLAB interface for the SPEX software package

    SPEX_Backslash  Exactly solve sparse linear systems with
                    default settings. This is the easiest
                    starting point for the SPEX software package.
                    SPEX_Backslash will automatically determine the
                    appropriate factorization algorithm for use in
                    solving your problem A x = b

    SPEX_Symmetric  Sparse integer-preserving SPEX_Symmetric
                    factorization for exactly solving SPD
                    linear systems, or symmetric indefinite
                    systems (with nonsingular minors)

    SPEX_LU         Sparse left-looking integer-preserving
                    LU factorization for exactly solve
                    sparse linear systems.

    SPEX_QR         Sparse integer-preserving QR factorization (FUTURE)

    SPEX_LU_ColRep  Sparse column replacement for SPEX LU factorization

    SPEX_Symmetric_Rank1: rank 1 updates for the SPEX Cholesky and LDL
                    factorizations

    SPEX_Utilities  Utility functions for all SPEX components

    SPEX_Update_Utilities   Utility functions for update/downdate
                    methods

    Makefile        compiles SPEX and its dependencies

    Python          SPEX python interface

Dependencies:

    AMD                 approximate minimum degree ordering

    COLAMD              column approximate minimum degree ordering

    SuiteSparse_config  configuration for all of SuiteSparse

    GNU GMP             GNU Multiple Precision Arithmetic Library
                        for big integer operations.  v6.1.2 or later
                        is required.

    GNU MPFR            GNU Multiple Precision Floating-Point Reliable
                        Library for arbitrary precision floating point
                        operations. v4.0.2 or later is required.

Compilation options:

* `SPEX_USE_PYTHON`:

  If `ON`, build Python interface for SPEX.
  If `OFF`: do not build the SPEX Python interface.
  Default: `SUITESPARSE_USE_PYTHON` (defaults to ON).

* `SPEX_USE_OPENMP`:

  If `ON`, OpenMP is used in SPEX if it is available.
  Default: `SUITESPARSE_USE_OPENMP` (defaults to ON).

To compile SPEX and its dependencies, just type "make" in this folder.
This will also run a few short demos.
To install the package system-wide, do "sudo make install"

Authors (alphabetical order):

    Jinhao Chen
    Timothy A. Davis
    Christopher Lourenco
    Lorena Mejia-Domenzain
    Erick Moreno-Centeno

