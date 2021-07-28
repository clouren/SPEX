function [zmax,PHIiter,PHIIiter,xbasic,ibasic]=linprog(A,b,c)
%
%LINPROG uses the two phase simplex method to solve the linear 
%program maximize cx subject to the constraints Ax = b and x >= 0 , 
%where A is m x n , rank(A) = m , and b >= 0 .    The output vector 
%is [zmax,PHIiter,PHIIiter,xbasic,ibasic], where zmax is the maximal 
%objective value, PHIiter and PHIIiter are the phase I and phase II 
%iteration counts, respectively, where xbasic is the vector of basic 
%x variables at optimality, and where ibasic is the indices of the 
%optimal basis columns in A (and hence the indices for the entries 
%in xbasic).  LINPROG detects infeasibility and unboundedness, and
%provides appropriate output messages in such cases.  LINPROG also
%contains a heuristic check for cycling, terminating the algorithm
%when m Phase II iterations occur without a change in the objective
%value.  See also PHASEI and PHASEII.
%
%Written for MATLAB version 5.0
%
%Written by Jeff Stuart, Department of Mathematics,
%University of Southern Mississippi, Hattiesburg, MS 39406.
%December, 1993.  Revised, October, 1997.
%jeffrey.stuart@usm.edu
%
[m,n]=size(A);
if max(size(b) ~=[m 1])
   disp('The dimensions of b do not match the dimensions of A.');
   return;
end   
if min(b) < 0
   disp('The RHS vector b must be nonnegative.');
   return;
end   
if max(size(c) ~=[1 n])
   disp('The dimensions of c do not match the dimensions of A.');
   return;
end
if rank(A) ~=m
   disp('A does not have full row rank.');
   return;
end
f=flops;
t=cputime;
PHIiter=0;
PHIIiter=0;
tol=0.0000000001;
xbasic=zeros(1,n);
[wmax,ibasic,PHIiter]=phasei(A,b);
if wmax < -tol
     f=flops -f;     
     t=cputime -t;
     disp('The original LP is infeasible.  Infeasibility was');
     disp('detected during Phase I.  The total number of phase');
     disp('one iterations performed was: '), disp(PHIiter);
     disp('The total flops required was: '),disp(f);
     disp('The required cpu time was: '),disp(t);
else
     disp('Phase I completed.  Original LP is feasible.');
     disp('The total number of Phase I iterations was: '),disp(PHIiter);
     disp('Starting Phase II.');
     [zmax,xbasic,ibasic,ienter,PHIIiter,PCOL,OPTEST,CYCTEST]=phaseii(A,b,c,ibasic);
     xbasic=xbasic';
     f=flops - f;
     t=cputime -t;
     if CYCTEST==1
        return;
     end
     if OPTEST == 0
          disp('The orginal LP is unbounded. An unbounded ray was');
          disp('detected during Phase II.  The output objective');
          disp('value is for the last basic solution found.');
          disp('The number of Phase II iterations was: '),disp(PHIIiter);
          disp('Last objective value is '),disp(zmax);
          disp('The last basic solution, xbasic is '),disp(xbasic);
          disp('The column indices for the last basis: '),disp(ibasic);
          disp('The index of the unbounded entering variable: '),disp(ienter);
          disp('The unbounded ray column is: '),disp(PCOL);
          disp('The total number of flops required was: '),disp(f);
          disp('The required cpu time was: '),disp(t);
     else
          disp('The original LP has an optimal solution.');
          disp('The number of Phase II iterations was: '),disp(PHIIiter);
          disp('The optimal objective value is '),disp(zmax);
          disp('The indices for the basic columns: '),disp(ibasic);
          disp('The optimal, basic solution is '),disp(xbasic);
          disp('The total number of flops required was: '),disp(f);
          disp('The required cpu time was: '),disp(t);
     end
end
