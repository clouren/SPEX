list = dir('LPnetlib/*.txt');
close all

%for i = 1:length(list)
for i = 32
    list(i).name
    fRead = fopen(strcat('LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);

    x_axis = 1:size(A,2);

    figure(i);
    title('number of nonzero');
    hold on
    plot(x_axis,A(4,:),'r',x_axis,A(6,:),'r--o');
    plot(x_axis,A(9,:),'b',x_axis,A(11,:),'b--o');
    plot(x_axis,A(14,:),'g',x_axis,A(16,:),'g--o');
    legend('L_{LLU}','U_{LLU}','L_{lb}','U_{lb}','L_{LUU}','U_{LUU}');

    figure(i+1);
    title('number of digit');
    hold on
    plot(x_axis,A(5,:),'r',x_axis,A(7,:),'r--o');
    plot(x_axis,A(10,:),'b',x_axis,A(12,:),'b--o');
    plot(x_axis,A(15,:),'g',x_axis,A(17,:),'g--o');
    legend('L_{LLU}','U_{LLU}','L_{lb}','U_{lb}','L_{LUU}','U_{LUU}');

    figure(i+2);
    title('factorization/updating time over updating time');
    hold on
    semilogy(x_axis,A(3,:)./A(13,:),'r');
    semilogy(x_axis,A(8,:)./A(13,:),'b');
    semilogy(x_axis,A(13,:)./A(13,:),'g');
    legend('t_{LLU}','t_{lb}','t_{LUU}');

    figure(i+3);
    title('solving/searching time');
    hold on
    plot(x_axis,A(1,:));
    plot(x_axis,A(2,:),'--');
    legend('solving 3 linear equations','search entering and existing column');
end
