function [z,xbasic,ibasic,ienter,iter,PCOL,OPTEST,CYCTEST]=phaseii(A,b,c,ibasic);
%PHASEII performs phase II of the simplex method starting with the basic
%columns specified by the vector ibasic.
%
%See also PHASEI and LINPROG.
%Written for Matlab version 5.0.
%
%Written by Jeff Stuart, Department of Mathematics, University of Southern Mississippi,
%Hattiesburg, MS 39406. October, 1993. Revised October, 1997.
%jeffrey.stuart@usm.edu
[m,n]=size(A);
PCOL=[];
ienter=[];
iter=0;
cycle=0;
CYCTEST=0;
X=zeros(1,n);
J=X;
J(ibasic)=ones(1,m);
K=[1:n];
inon=K(~J);
B=A(:,ibasic);
xbasic=B\b;
z=c(ibasic)*xbasic;
if m<n;
   X(ibasic)=xbasic;
   Cred=c(inon) - (c(ibasic)/B)*A(:,inon);
   OPTEST=1;
   loop =1;
   while loop ==1;
      if max(Cred) > 0; 
         iter=iter + 1;
         [Maxcost,j]=max(Cred);
         ienter=inon(j);
         PCOL=B\A(:,ienter);
         if PCOL <= 0 , OPTEST = 0;
            loop = 0;
         else
            J(ienter)=1;
            TESTROWS=find(PCOL > 0);
            TESTCOL=PCOL(TESTROWS);
            [minrat,j]=min(xbasic(TESTROWS)./TESTCOL);
            if minrat <=0, cycle = cycle+1;
               if cycle > m;
                  disp('Algorithm terminated due to excessive cycling.')
                  disp('Restart algorithm from phase II using a perturbed')
                  disp(' RHS vector b and the current basis.')
                  disp(ibasic)
                  CYCTEST=1;
                  break
               end
            else
               cycle = 0;
            end
            iexit=ibasic(TESTROWS(j));
            J(iexit)=0;
            xbasic=xbasic - minrat*PCOL;
            X(ibasic)=xbasic;
            X(ienter)=minrat;
            X(iexit)=0;
            z=z + Maxcost*minrat;
            J=logical(J);
            ibasic=K(J);
            inon=K(~J);
            B=A(:,ibasic);
            xbasic=X(ibasic)';
            Cred=c(inon) - (c(ibasic)/B)*A(:,inon);
         end
      else
        loop = 0;
      end
   end
end
