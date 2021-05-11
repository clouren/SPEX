function [wmax,ibasic,PHIiter]=phasei(A,b)
%PHASEI performs Phase I of the simplex method on the constraints 
%Ax = b and x >= 0 (where A is m x n, rank(A) = m , and b >= 0)
%to determine whether there is a feasible point.  The function 
%output is wmax, the artificial objective value; ibasic, the 
%indices of the basic variables at optimality; and PHIiter, the 
%number of Phase I iterations performed.  If  wmax < 0 ,then the 
%original LP is infeasible.  If wmax = 0 , the original LP is 
%feasible, and ibasic is the index set for a feasible basis.
%To allow for round-off error, the tests are  wmax < -tol  for
%infeasibility, and  wmax >= -tol  for feasibility, where "tol"
% is a preset tolerance (see the initialization value in the fourth
%line below).  At the expense of additional computation, an adaptive
%choice for "tol" based on A and b could be selected.
%
%See also PHASEII and LINPROG.
%Written for MATLAB version 5.0 .
%Written by Jeff Stuart, Department of Mathematics, University of 
%Southern Mississippi, Hattiesburg, MS 39406.  October, 1997.
%jeff.stuart@usm.edu
[m,n]=size(A);
A=[A,eye(m)];
PHIiter=0;
tol=0.0000001;
ztol=0.0000001;
X=zeros(1,n+m);
J=[zeros(1,n),ones(1,m)];
c=-J;
K=[1:n+m];
J=logical(J);
ibasic=K(J);
inon=K(~J);
B=eye(m);
xbasic=b;
w=-sum(xbasic);
X(ibasic)=b;
Cred= ones(1,m)*A(:,inon);
loop =1;
while loop ==1;
  if max(Cred) > ztol ;
    PHIiter=PHIiter + 1;
    [Maxcost,j]=max(Cred);
    ienter=inon(j);
    PCOL=B\A(:,ienter);
    J(ienter)=1;
    TESTROWS=find(PCOL > ztol);
    TESTCOL=PCOL(TESTROWS);
    [minrat,j]=min(xbasic(TESTROWS)./TESTCOL);
    iexit=ibasic(TESTROWS(j));
    J(iexit)=0;
    if minrat > 0;
       xbasic=xbasic - minrat*PCOL;
    end
    X(ibasic)=xbasic;
    X(ienter)=minrat;
    X(iexit)=0;
    w=w + Maxcost*minrat;
    ibasic=K(J);
    inon=K(~J);
    B=A(:,ibasic);
    xbasic=X(ibasic)';
    Cred=c(inon) - (c(ibasic)/B)*A(:,inon);
  elseif Cred <= ztol;
    loop = 0;
  end
end
wmax=-sum(X(n+1:n+m));
if wmax >= -tol;
   X=X(1:n);
   last=ibasic(m);
   K=K(1:last);
   ibasic=K(J(1:last));
   while last > n;
      J=J(1:last);
      K=[1:last];
      inon=K(~J);
      B=A(:,ibasic);
      inon=inon(inon <= n);
      j=find(([zeros(1,m-1),1]/B)*A(:,inon));
      ienter=inon(j(1));
      J(ienter)=1;
      J(last)=0;
      ibasic=K(J);
      last=ibasic(m);
      PHIiter=PHIiter+1;
   end
end
