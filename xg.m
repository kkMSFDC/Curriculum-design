function R1 = xg(m1,m2)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
N=length(m1);
t1=m1;
t2=zeros(1,N-1);
R1=[];
for i=1:N+1
    D=sum(1*(t1~=m2));
    A=N-D;
    R1(i)=(A-D)/N;
    for j=1:N-1
        t2(j)=t1(j);
    end
    flag=t1(N);
    for j=1:N-1
        t1(j+1)=t2(j);
    end
    t1(1)=flag;
end
end
    

