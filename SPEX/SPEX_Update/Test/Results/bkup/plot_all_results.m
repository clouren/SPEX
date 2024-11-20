list = dir('LPnetlib/*.txt');
close all

max_ratio = 0;
for i = 1:length(list)
    list(i).name
    fRead = fopen(strcat('LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);

    t1(i) = sum(A(3,2:end));
    t2(i) = sum(A(8,2:end));
    t3(i) = sum(A(13,2:end));
    for j = 1:size(A,2)
        if A(13,j)~=0 && max_ratio < A(3,j)/A(13,j)
            max_ratio = A(3,j)/A(13,j)
            A(3,j)
            A(13,j)
        end
    end
    max_ratio
end
figure(1);
semilogy(1:length(list),t1./t3,'r-o');
title('factorization/updating time over updating time');
hold on
semilogy(1:length(list),t2./t3,'b-o');
semilogy(1:length(list),t3./t3,'g-o');
legend('t_{LLU}','t_{lb}','t_{LUU}');

figure(2);
title('factorization/updating time over updating time');
hold on
plot(1:length(list),t1./t3,'r');
plot(1:length(list),t2./t3,'b');
plot(1:length(list),t3./t3,'g');
legend('t_{LLU}','t_{lb}','t_{LUU}');
