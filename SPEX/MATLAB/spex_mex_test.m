function spex_mex_test
% SPEX_MEX_TEST: run a set of tests for SPEX matlab interface
%
% Usage:  spex_mex_test
%
% See also spex_backslash, spex_lu_backslash, spex_cholesky_backslash,
%   spex_ldl_backslash, spex_mex_install, spex_mex_demo.

% SPEX: (c) 2022-2024, Christopher Lourenco, Jinhao Chen,
% Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
% All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

maxerr = 0 ;
rng ('default') ;

fprintf ('Testing SPEX Left LU: ') ;

% First, check if we can use a real life sparse matrix via ssget
if (exist ('ssget') ~= 0)
    fprintf ('. (please wait) ') ;
    % 159 is a square SPD matrix
    prob = ssget(159);
    A = prob.A;
    [m n] = size(A);
    b = rand(m, 1);
    fprintf ('.') ;
    x = spex_lu_backslash(A,b);
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;

    % now convert to an integer problem (x will not be integer)
    A = floor (2^20 * A) ;
    b = floor (2^20 * b) ;
    fprintf ('.') ;
    x = spex_lu_backslash (A, b) ;
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;
    fprintf ('.') ;
end

orderings = { 'default', 'none', 'colamd', 'amd' } ;
pivotings = { 'smallest', 'diagonal', 'first', ...
    'tol smallest', 'tol largest', 'largest' } ;

for n = [1 10 100]
    for density = [0.001 0.05 0.5]

        % construct a well-conditioned problem to solve
        A = sprand(n,n,density);
        A = A+A' + n * speye (n) ;
        b = rand(n,1);

        for korder = 1:length (orderings)
            for kpiv = 1:length (pivotings)
                for tol = [0.1 0.5]

                    clear option
                    option.order = orderings {korder} ;
                    option.pivot = pivotings {kpiv} ;
                    option.tol   = tol ;

                    fprintf ('.') ;
                    x = spex_lu_backslash(A,b, option);
                    x2 = A\b;
                    err = norm(x-x2)/norm(x);
                    maxerr = max (maxerr, err) ;

                    % now convert to an integer problem (x will not be integer)
                    A = floor (2^20 * A) ;
                    b = floor (2^20 * b) ;
                    x = spex_lu_backslash(A,b, option);
                    x2 = A\b;
                    err = norm(x-x2)/norm(x);
                    maxerr = max (maxerr, err) ;
                end
            end
        end
    end
end

fprintf ('\nmaxerr: %g\n', maxerr) ;

if (maxerr < 1e-6)
    fprintf('\nSPEX LU tests successful\n')
else
    error ('\nTesting failure!  error too high\n')
end

fprintf ('Testing SPEX Cholesky and LDL: ') ;

% First, check if we can use a real life sparse matrix via ssget
if (exist ('ssget') ~= 0)
    fprintf ('. (please wait) ') ;
    % 2 is a square SPD matrix
    prob = ssget(2);
    A = prob.A;
    [m n] = size(A);
    b = rand(m, 1);
    fprintf ('.') ;
    x = spex_cholesky_backslash(A,b);
    x1 = spex_ldl_backslash(A,b);
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;
    err = norm(x-x1)/norm(x);
    maxerr = max (maxerr, err) ;

    % now convert to an integer problem (x will not be integer)
    A = floor (2^20 * A) ;
    b = floor (2^20 * b) ;
    fprintf ('.') ;
    x = spex_cholesky_backslash (A, b) ;
    x1 = spex_ldl_backslash (A, b) ;
    x2 = A\b;
    err = norm(x-x1)/norm(x);
    maxerr = max (maxerr, err) ;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;

    % test a negative definite matrix with LDL
    x3 = spex_ldl_backslash (-A, b) ;
    x4 = (-A)\b;
    err = norm(x3-x4)/norm(x);
    maxerr = max (maxerr, err) ;
    fprintf ('.') ;

end

orderings = { 'none', 'colamd', 'amd' } ;
%pivotings = { 'smallest', 'diagonal', 'first', ...
%    'tol smallest', 'tol largest', 'largest' } ;

for n = [1 10 100]
    for density = [0.001 0.05 0.5]

        % construct a well-conditioned problem to solve
        A = sprand(n,n,density);
        A = A+A' + n * speye (n) ;
        b = rand(n,1);

        for korder = 1:length (orderings)

            clear option
            option.order = orderings {korder} ;

            fprintf ('.') ;
            x = spex_cholesky_backslash(A,b, option);
            x1 = spex_ldl_backslash(A,b, option);
            x2 = A\b;
            err = norm(x-x1)/norm(x);
            maxerr = max (maxerr, err) ;
            err = norm(x-x2)/norm(x);
            maxerr = max (maxerr, err) ;

            % now convert to an integer problem (x will not be integer)
            A = floor (2^20 * A) ;
            b = floor (2^20 * b) ;
            x = spex_cholesky_backslash(A,b, option);
            x1 = spex_ldl_backslash(A,b, option);
            x2 = A\b;
            err = norm(x-x1)/norm(x);
            maxerr = max (maxerr, err) ;
            err = norm(x-x2)/norm(x);
            maxerr = max (maxerr, err) ;
        end
    end
end

fprintf ('\nmaxerr: %g\n', maxerr) ;

if (maxerr < 1e-6)
    fprintf('\nSPEX Cholesky and LDL tests successful\n')
else
    error ('spex_cholesky_backslash:test', '\nTesting failure!\n')
end


fprintf ('Testing SPEX Backslash: ') ;
for n = [1 10 100]
    for density = [0.001 0.05 0.5]

        % construct a well-conditioned problem to solve
        A = sprand(n,n,density);
        A = A + n*speye(n);
        A2 = A+A' + n * speye (n) ;
        b = rand(n,1);

        for korder = 1:length (orderings)

            clear option
            option.order = orderings {korder} ;

            fprintf ('.') ;
            x = spex_cholesky_backslash(A2,b, option);
            x1 = spex_ldl_backslash(A2,b, option);
            x2 = A2\b;
            err = norm(x-x1)/norm(x);
            maxerr = max (maxerr, err) ;
            err = norm(x-x2)/norm(x);
            maxerr = max (maxerr, err) ;

            fprintf ('.')
            x = spex_lu_backslash(A, b, option);
            x2 = A\b;
            err = norm(x-x2)/norm(x);
            maxerr = max( maxerr, err);

            % now convert to an integer problem (x will not be integer)
            A2 = floor (2^20 * A2) ;
            b = floor (2^20 * b) ;
            x = spex_cholesky_backslash(A2,b, option);
            x1 = spex_ldl_backslash(A2,b, option);
            x2 = A2\b;
            err = norm(x-x1)/norm(x);
            maxerr = max (maxerr, err) ;
            err = norm(x-x2)/norm(x);
            maxerr = max (maxerr, err) ;

            % now convert to an integer problem (x will not be integer)
            A = floor (2^20 * A) ;
            x = spex_lu_backslash(A,b, option);
            x2 = A\b;
            err = norm(x-x2)/norm(x);
            maxerr = max (maxerr, err) ;
        end
    end
end

if (maxerr < 1e-6)
    fprintf('\nSPEX Backslash tests\n')
else
    error ('SPEX_backslash:test', '\nTesting failure!  error too high. Please reinstall\n')
end

fprintf("\nAll testing complete, ready to go!\n");
