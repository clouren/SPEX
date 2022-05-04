close all
clc
clear all
list = dir('LPnetlib_CholUpdate/*.txt');

max_ratio = 0;
min_ratio = Inf;
k = length(list);
nodata = 0;
too_simple = 0;
too_simple1 = 0;
too_simple2 = 0;
for i = 1:length(list)
    i
    list(i).name
    fRead = fopen(strcat('LPnetlib_CholUpdate/',list(i).name), 'r');
    A = fscanf(fRead, '%f %d %d %f %d %d %f %d %d %f %d %d %f %d %d %f %d %d %f %d %d',[21, 1]);
    fclose(fRead);
    
    t_fact(i) = 0;
    t_update(i) = 0;
    if (length(A) ~= 0)
        t_fact(i) = A(1)*2+A(7)+A(16);
        t_update(i) = A(4)+A(10)+A(13)+A(19);
    else
        nodata = nodata+1
    end
%     figure(2);
%     loglog(t_fact,t_update,'b*');
%     hold on;
%     if i>1
%         loglog(t_fact(1:i-1),t_update(1:i-1),'bo');
%     end
%     loglog(t_fact(i),t_update(i),'r*');
%     if i == 69 || i ==70
%         figure(i);
%         loglog(t_fact,t_update,'bo');
%         hold on;
%         loglog(t_fact(i),t_update(i),'ro');
%     end

    if t_update(i) == 0
        too_simple1 = too_simple1+1
    end
    if t_fact(i) == 0
        too_simple2 = too_simple2+1
    end
    if t_update(i) == 0 || t_fact(i) == 0
        too_simple = too_simple+1;
    end
    i-too_simple
    
    if t_update(i)~=0 && t_fact(i)/t_update(i) > max_ratio
        max_ratio =t_fact(i)/t_update(i)
    end
    if t_fact(i)~=0 && t_fact(i)/t_update(i) < min_ratio
        min_ratio =t_fact(i)/t_update(i)
    end
    
end
figure(1);
%title('factorization/updating time over updating time');
xticks(0: 10: 80);
xlabel("case index");
ylabel("time ratio");
hold on
semilogy(1:k,t_fact./t_update,'b-o');
semilogy(1:k,t_update./t_update,'g-o');
legend('t_{DLU}','t_{LUU}');

figure(3);
%title('factorization/updating time over updating time');
loglog(t_fact,t_fact);
hold on
loglog(t_fact,t_fact/10);
loglog(t_fact,t_fact/100);
loglog(t_fact,t_update,'bo');
xlabel("Time for direct Cholesky factorization (sec)");
ylabel("Time for Cholesky update/downdate (sec)");
legend('1x','1/10x','1/100x');