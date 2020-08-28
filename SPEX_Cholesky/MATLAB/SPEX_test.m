function IP_test
%IP_test: run a set of tests for IP_Chol
%
% Usage:  IP_test

maxerr = 0 ;

fprintf ('Testing IP_Chol: ') ;

% First, check if we can use a real life sparse matrix via ssget
if (exist ('ssget') ~= 0)
    fprintf ('.') ;
    % 159 is a square SPD matrix
    prob = ssget(159);
    A = prob.A;
    [m n] = size(A);
    b = rand(m, 1);
    maxerr = max (maxerr, IP_run_all_tests (A, b)) ;
end

for n = [1 10 100]
    for density = [0.001 0.05 0.5 1]
        fprintf ('.') ;
        A = sprand(n,n,density);
        A = A'*A;
        % Want a numerically stable A
        if (condest(A) > 1e6)
            A = A + speye (n) ;
        end
        b = rand(n,1);
        maxerr = max (maxerr, IP_run_all_tests (A, b)) ;
    end
end

fprintf ('\nmaxerr: %g\n', maxerr) ;

if (maxerr < 1e-10)
    fprintf('\nTesting complete, installation successful\n')
else
    error ('IP_Chol:test', '\nTesting failure!  error too high\n')
end
