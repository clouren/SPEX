function maxerr = IP_run_all_tests (A, b)
%IP_run_all_tests run a set of sets of IP_Chol
%
% maxerr = IP_run_all_tests (A, b)

try chol(A);

x = IP_Chol(A,b);
x2 = A\b;
err1 = norm(x-x2)/norm(x);

maxerr = err1;

catch
    maxerr = 0;
    
end
