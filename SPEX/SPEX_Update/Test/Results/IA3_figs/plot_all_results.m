close all
clc
clear all
list = dir('LPnetlib/*.txt');

max_ratio = 0;
k = 0;
for i = 1:length(list)
    %list(i).name
    fRead = fopen(strcat('LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);
    
    if sum(A(13,2:end))~=0
        k=k+1;
        t1(k) = sum(A(3,2:end));
        t2(k) = sum(A(8,2:end));
        t3(k) = sum(A(13,2:end));
    else
        i
    end
    for j = 1:size(A,2)
        if A(13,j)~=0 && max_ratio < A(3,j)/A(13,j)
            max_ratio = A(3,j)/A(13,j)
            A(3,j)
            A(13,j)
        end
    end
    %max_ratio
end
figure(1);
semilogy(1:k,t1./t3,'r-o');
%title('factorization/updating time over updating time');
xticks(0: 10: 80);
xlabel("case index");
ylabel("time ratio");
hold on
semilogy(1:k,t2./t3,'b-o');
semilogy(1:k,t3./t3,'g-o');
legend('t_{DLU}','t_{lb}','t_{LUU}');

figure(3);
%title('factorization/updating time over updating time');
loglog(t1,t1);
hold on
loglog(t1,t1/10);
loglog(t1,t1/100);
loglog(t1,t3,'bo');
xlabel("Time for direct LU factorization (sec)");
ylabel("Time for LU update (sec)");
legend('1x','1/10x','1/100x');