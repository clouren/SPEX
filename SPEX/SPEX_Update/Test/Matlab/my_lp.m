clc
close all
clear all
Prob=ssget(594)
A=Prob.A;
c=Prob.aux.c;
b=Prob.b;
lb=Prob.aux.lo;
ub=Prob.aux.hi;

fval_add=0;
%Check if any variables have equal upper and lower bounds. 
%If so, check for feasibility, and then fix and remove the variables.
SAME_BOUND =(lb==ub);
if sum(SAME_BOUND)>0
    b=b-A(:,SAME_BOUND)*SAME_BOUND*lb;
    fval_add=fval_add+c'*SAME_BOUND*lb;
    c(SAME_BOUND)=[];
    A(:,SAME_BOUND)=[];
    lb(SAME_BOUND)=[];
    ub(SAME_BOUND)=[];
    fprintf("1:delete columns of A:\n");
    find(SAME_BOUND)
end

%Check if any linear inequality constraint involves only one variable.
%If so, check for feasibility, and then change the linear constraint to a bound.
%skip since no inequality

%Check if any linear equality constraint involves only one variable.
%If so, check for feasibility, and then fix and remove the variable.
i = 1;
while i <= size(A,1)
    if nnz(A(i,:))==1
        [~,j,s]=find(A(i,:));
        xj = b(i)/s;
        if (xj>=lb(j) && xj<=ub(j))
            fval_add=fval_add+c(j)*xj;
            b=b-A(:,j)*xj;
            c(j)=[];
            A(i,:)=[];
            A(:,j)=[];
            b(i)=[];
            lb(j)=[];
            ub(j)=[];
            fprintf("2:delete col of A :%d\n",j);
            fprintf("2:delete row of A :%d\n",i);
        else
            error("infeasible condition\n");
        end
    else
        i = i+1;
    end
end

%Check if any linear constraint matrix has zero rows. If so, check for feasibility, and then delete the rows.
i = 1;
while i<=size(A,1)
    if nnz(A(i,:)) == 0
        if b(i) == 0
            b(i) = [];
            A(i,:)=[];
            fprintf("3:delete zero row of A :%d\n",i);
        else
            error("infeasible condition\n");
        end
    else
        i = i+1;
    end
end

%Determine if the bounds and linear constraints are consistent.

%Check if any variables appear only as linear terms in the objective function and do not appear in any linear constraint.
%If so, check for feasibility and boundedness, and then fix the variables at their appropriate bounds.
j = 1;
while j<=size(A,2)
    if nnz(A(:,j)) == 0
        if c(j)>0
            fval_add=fval_add+c(j)*lb(j);
        elseif c(j) < 0
            fval_add=fval_add+c(j)*ub(j);
        end
        c(j)=[];
        A(:,j)=[];
        lb(j)=[];
        ub(j)=[];
        fprintf("4:delete col of A :%d\n",j);
    else
        j = j+1;
    end
end

%Change any linear inequality constraints to linear equality constraints by adding slack variables.

%check if a column is simply multiple of another one 
for i =size(A,2)-1:-1:1
    if (i==305)
        fprintf("here");
    end
    A_coli_1stnz=find(A(:,i),1);
    J=i+find(A(A_coli_1stnz,i+1:size(A,2)));
    for j_index=length(J):-1:1
        j=J(j_index);
        if nnz(A(:,i)) == nnz(A(:,j))
            Ai_nz_pattern =find(A(:,i));
            Aj_nz_pattern =find(A(:,j));
            if isempty(find(Ai_nz_pattern~=Aj_nz_pattern,1))
                %same pattern
                multipliers = A(:,j)./A(:,i);
                if nnz(A(:,i)) ==1 || isempty(find(multipliers~=multipliers(1),1))
                    A(:,i)
                    A(:,j)
                    A(:,j)=[];
                    c(j)=[];
                    if multipliers(1) > 0
                        lb(i)=lb(i)+multipliers(1)*lb(j);
                        ub(i)=ub(i)+multipliers(1)*ub(j);
                    else
                        lb(i)=lb(i)+multipliers(1)*ub(j);
                        ub(i)=ub(i)+multipliers(1)*lb(j);
                    end
                    lb(j)=[];
                    ub(j)=[];
                    fprintf("5:delete col of A :%d, which is multiple of col %d of A\n",j,i);
                end
            end
        end
    end
end

for i=1:size(A,2)
    col_nnz(i)= nnz(A(:,i));
end
selected=(col_nnz==1);
for j=1:size(A,1)
    row_nnz(j)=nnz(A(j,selected));
end
used_row=(row_nnz==1);
%add column of all -1s if min(b) < 0
A_aux=A;
if min(b) < 0
    A_aux=[A, -1*ones(size(A,1),1)];
    [~,minb_index]=min(b);
    if nnz(A(minb_index,selected))==1
        for j=1:size(A,2)
            if col_nnz(j)==1 && find(A(:,j),1)==minb_index
                selected(j)=0;
            end
        end
    end
    selected(size(A_aux,2))=1;
    used_row(minb_index)=1;
end

for j=1:size(A,1)
    row_nnz(j)=nnz(A(j,selected==0));
end

[~,index]=min(row_nnz(used_row==0));
index = find(used_row==0,index);
i=index(end);
while i<=size(A,1)
    J=find(A_aux(i,:));
    for j=J
        if selected(j)==0
            selected(j)=1;
            used_row(i)=1;
            col_nnz_pattern=find(A_aux(:,j));
            row_nnz(col_nnz_pattern)=row_nnz(col_nnz_pattern)-1;
            if sum(used_row)< size(A,1)
                [~,index]=min(row_nnz(used_row==0));
                index = find(used_row==0,index);
                i=index(end);
            else
                i=size(A,1)+1;
            end
            break;
        end
    end
end

B = A_aux(:,selected==1);
[sum(selected) sum(used_row) sprank(B)]
[L1,U1]=lu(B);
find(diag(U1)==0)
x(selected==1)=B\b;
norm(A_aux*x'-b)
[p,q,~,~,~,~] = dmperm(A); C=A(p,q);
figure(2)
spy(C)
selected=zeros(1,size(A,2));
selected(q(size(A,2)-size(A,1)+1:end))=1;
%selected(q(1:size(A,1)))=1;
B = A(:,selected==1);
sprank(B)
[L1,U1]=lu(B);
find(diag(U1)==0)

x=zeros(size(A,2),1);

x(selected==1)=B\b;
norm(A*x-b)

[L,U]=lu(A);
[p,q,r,s,cc,rr] = dmperm(U); C=U(p,q);
figure(3)
spy(U)
figure(4)
spy(C)

selected=zeros(1,size(A,2));
selected(q(size(A,2)-size(A,1)+1:end))=1;
%selected(q(1:size(A,1)))=1;
B = A(:,selected==1);
sprank(B)
[L1,U1]=lu(B);
find(diag(U1)==0)

x=zeros(size(A,2),1);

x(selected==1)=B\b;
norm(A*x-b)