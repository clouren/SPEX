
Prob=ssget(597);
A=Prob.A;
c=Prob.aux.c;
b=Prob.b;
lb=Prob.aux.lo;
ub=Prob.aux.hi;
%options = optimoptions('linprog','Algorithm','dual-simplex','Display','iter')
%[x,fval,exitflag,output,lambda] =LinProg(c,[],[],A,b,lb,ub,options); c'*x

basis = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 39 15 16 17 18 21 22 32 33 37 38 49 50]+1;
%basis = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 39 15 16 17 18 21 22 23 32 33 37 38 49]+1;
x_sol=(A(:,basis)\b)
c(basis)'*x_sol;
used = zeros(1,51);
used(basis) = 1;
c_updated=(A(:,used==1)'\c(used==1))'*A(:,used==0)-c(used==0)'
unused_list=find(used==0)
[~,max_ind]=max(c_updated);
y_sol=A(:,basis)\A(:,unused_list(max_ind))
x_sol(y_sol~=0)./y_sol(y_sol~=0)

A_updated=A(:,basis)\full(A)
c_new=c';
for p= 1:27
    i = basis(p);
    if c_new(i) ~= 0
        c_new = c_new-c_new(i)*A_updated(p,:)
    end
end


