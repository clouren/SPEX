
Prob=ssget(594)
A=Prob.A;
c=Prob.aux.c;
b=Prob.b;
lb=Prob.aux.lo;
ub=Prob.aux.hi;
options = optimoptions('linprog','Algorithm','dual-simplex','Display','iter')
[x,fval,exitflag,output,lambda] =LinProg(c,[],[],A,b,lb,ub,options); c'*x