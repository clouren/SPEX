%% a demo of the SPEX QR MATLAB interface
% SPEX is a package for solving sparse linear systems of
% equations with a roundoff-free integer-preserving method.
% The result is always exact, unless the matrix A is perfectly
% singular.
%
% See also vpa, spex_backslash, spex_lu_backslash,
%   spex_cholesky_backslash, spex_ldl_backslash,
%   spex_mex_install, spex_mex_test.
%
% Copyright (c) 2022-2026, Christopher Lourenco, Jinhao Chen,
% Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
% All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

%#ok<*NOPTS>
%#ok<*NASGU>

% SPEX QR can exactly solve Ax = b if A is square and full rank
% If A is rectangular with m > n, SPEX QR returns an exact least
% squares solution if A has full column rank or an exact basic solution
% if A is rank deficient.
% If A is rectangular with n > m, SPEX QR returns the exact minimum norm
% solution if A has full row rank and an exact basic solution if A is rank
% deficient.
% Finally, the spex_rank function can return the exact rank of a matrix.
% This function utilizes LU factorization if A is square or QR
% factorization if A is rectangular. In this function, we demo each of
% these functionalities.

%% SPEX QR exact solutions of square full rank matrices
% If A is square and full rank, SPEX QR returns the same solution as an LU
% or Cholesky solver. In this case, it is recommended to use LU or Cholesky
% instead; however, QR can return exact solutions.
A = rand(50);
b = ones(50,1);
x = spex_qr_backslash(A,b);
x2 = spex_lu_backslash(A,b);
% This should be zero because both methods return the same solution
err_spex_square = norm(x-x2);

%% SPEX QR solving rectangular systems with m > n
% If m is larger than n, we obtain an exact least squares solution if A is
% full column rank and a basic solution otherwise. We start with a few full
% rank systems. Note here backslash and spex should return the same least
% squares solution.
A = rand(30,20); b = ones(30,1);
x = spex_qr_backslash(A,b);
x2 = spex_cholesky_backslash(A'*A,A'*b);
x3 = A\b;
err_spex_ls = norm(x-x2);
err_matlab_ls = norm(x-x3);

% Now we will modify A to be rank deficient
A(:,[2 5 7]) = ones(30,3);
% This will return a basic solution
x = spex_qr_backslash(A,b);
x2 = A\b;
spex_rankdef_norm1 = norm(A*x-b);
matlab_rankdef_norm1 = norm(A*x2-b);
diff_norms = spex_rankdef_norm1 - matlab_rankdef_norm1;

% One more time with a different matrix
A = rand(50,25); A(:,3) = A(:,7); A(:,13) = A(:,17); A(:,23) = A(:,25);
b = ones(50,1);
x = spex_qr_backslash(A,b);
x2 = A\b;
spex_rankdef_norm2 = norm(A*x-b);
matlab_rankdef_norm2 = norm(A*x2-b);
diff_norms2 = spex_rankdef_norm2 - matlab_rankdef_norm2;

%% SPEX QR solving rectangular systems with n > m
% If n is larger than m, SPEX will return a minimum norm solution. This is
% different than the solution returned by matlab backslash, however we can
% compare with matlab's built in least squares min norm solution. As usual,
% we will start with full rank and then try rank deficient.
A = rand(20,30); b = ones(20,1);
x = spex_qr_backslash(A,b);
x2 = lsqminnorm(A,b);
err_matlab_lsqnorm = norm(x-x2);

% Now we try rank deficient matrices
A(3,:) = A(6,:); A(5,:) = A(8,:); A(10,:) = A(11,:); A(14,:) = A(20,:);
x = spex_qr_backslash(A,b);
x2 = lsqminnorm(A,b);
spex_axb_norm = norm(A*x-b);
matlab_axb_norm = norm(A*x2 - b);
err_matlab_lsqnormdeficient = norm(x-x2);


%% SPEX Rank computing the exact rank of different matrices

%% Print output
fprintf("\nSPEX QR with full rank square linear systems")
fprintf("\nSPEX QR vs SPEX LU")
err_spex_square

fprintf("\nSPEX QR with full rank m >= n linear systems")
fprintf("\nQR LS vs Chol normal equations")
err_spex_ls
fprintf("\nQR vs MATLAB least squares")
err_matlab_ls

fprintf("\nSPEX QR with rank deficient m >= n linear systems")
fprintf("\nQR norm Ax-b matrix 1")
spex_rankdef_norm1
fprintf("\nMATLAB norm Ax-b matrix 1")
matlab_rankdef_norm1

fprintf("\nQR norm Ax-b matrix 2")
spex_rankdef_norm2
fprintf("\nMATLAB norm Ax-b matrix 2")
matlab_rankdef_norm2

fprintf("\nSPEX QR with full rank n>m linear systems")
fprintf("\nQR vs matlab min norm")
err_matlab_lsqnorm

fprintf("\nSPEX QR with rank deficient n>m linear systems")
fprintf("\nQR accuracy norm Ax-b")
spex_axb_norm
fprintf("\nMATLAB accuracy norm Ax-b")
matlab_axb_norm
fprintf("QR vs matlab")
err_matlab_lsqnormdeficient