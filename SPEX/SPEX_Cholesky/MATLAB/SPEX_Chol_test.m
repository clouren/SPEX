function SPEX_Chol_test
%SPEX_Chol_test: run a set of tests for SPEX_Chol_backslash
%
% Usage:  SPEX_Left_LU_test
%
% See also SPEX_install, SPEX_backslash, SPEX_demo.

% SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
% Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
% All Rights Reserved.  See SPEX_Cholesky/License for the license.

maxerr = 0 ;
rng ('default') ;

fprintf ('Testing SPEX_Chol_backslash: ') ;

% First, check if we can use a real life sparse matrix via ssget
if (exist ('ssget') ~= 0)
    fprintf ('. (please wait) ') ;
    % 2 is a square SPD matrix
    prob = ssget(2);
    A = prob.A;
    [m n] = size(A);
    b = rand(m, 1);
    fprintf ('.') ;
    x = SPEX_Chol_backslash(A,b);
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;

    % now convert to an integer problem (x will not be integer)
    A = floor (2^20 * A) ;
    b = floor (2^20 * b) ;
    fprintf ('.') ;
    x = SPEX_Chol_backslash (A, b) ;
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;
    fprintf ('.') ;
end

orderings = { 'none', 'colamd', 'amd' } ;
%pivotings = { 'smallest', 'diagonal', 'first', ...
%    'tol smallest', 'tol largest', 'largest' } ;

for n = [1 10 100]
    for density = [0.001 0.05 0.5 1]

        % construct a well-conditioned problem to solve
        A = sprand(n,n,density);
        A = A+A' + n * speye (n) ;
        b = rand(n,1);

        for korder = 1:length (orderings)
            %for kpiv = 1:length (pivotings)
            %    for tol = [0.1 0.5]

                    clear option
                    option.order = orderings {korder} ;
%                    option.pivot = pivotings {kpiv} ;
%                    option.tol   = tol ;

                    fprintf ('.') ;
                    x = SPEX_Chol_backslash(A,b, option);
                    x2 = A\b;
                    err = norm(x-x2)/norm(x);
                    maxerr = max (maxerr, err) ;

                    % now convert to an integer problem (x will not be integer)
                    A = floor (2^20 * A) ;
                    b = floor (2^20 * b) ;
                    x = SPEX_Chol_backslash(A,b, option);
                    x2 = A\b;
                    err = norm(x-x2)/norm(x);
                    maxerr = max (maxerr, err) ;
                %end
            %end
        end
    end
end

fprintf ('\nmaxerr: %g\n', maxerr) ;

if (maxerr < 1e-6)
    fprintf('\nTesting complete, installation successful\n')
else
    error ('SPEX_Chol_backslash:test', '\nTesting failure!  error too high\n')
end
