clc
clear
list = dir('LPnetlib/*.txt');

k = 0;
for i = 1:length(list)
    fRead = fopen(strcat('LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);
    
    if sum(A(13,2:end))~=0
        k=k+1;
        t1(k) = sum(A(3,2:end));
        file(k) = i;
    end
end
[~,Index] = sort(t1);

for start=1:8
    start
    for i =start:8:k
        j=Index(i);
        slen=length(list(j).name);
        fprintf("./cholupdate %.*s > Outputs/CholUpdate/%.*s.out\n",slen-9,list(j).name,slen-4,list(j).name);
    end
    fprintf("\n\n\n");
end