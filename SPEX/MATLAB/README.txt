This is the first draft of a general SPEX_backslash function

SPEX_backslash_install installs each SPEX_* MATLAB interface and temporarily
adds them to the path. Then x = SPEX_backslash(A,b) will analyze A and provide
the appropriate function. Currently, if A appears to be symmetric, Cholesky is tried,
otherwise LU is done.
