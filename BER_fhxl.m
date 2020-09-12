%扩频系统
br=zeros(9,11);
NN=100000;
N=1;
v=1;
wco=0.5;
wc=0.8;
p=4;
Fs=100;
ts=1/Fs;

%产生7位m序列
m=m_se([0 1 1]);
m1=1-2*m;

%产生16位Walsh序列
wals=walsh_se(16);
wals1=wals(8,:);

%产生112位复合序列
CM=kron(wals1,m1);

gsr=10^(15/10);    %能量比值
sgma=sqrt(2*N*gsr);
t=ts:ts:112;
%产生单频干扰

%snr=0;
for f0=10:1:10
    yg=sgma*sin(2*pi*f0*t);
for snr=0:1:10
    BER1=0;
    for i=1:NN
        Fs=100;        %采样频率
        Fc=10;         %载波频率
        ts=1/Fs;       %采样间隔
        fd=1;          %码元速率
        N=1;         %码元个数

        %产生信息码元
        code2=randi([0,1],[1,N]);
        code3=2.*code2-1;

       
        %[m2,nn]=pulse(m1,Fs,fd);         %脉冲信号转换为矩形波

         %扩频
         L=kron(code3,CM);
         
         %平滑
         [B,A]=butter(p,wc);
         L=filtfilt(B,A,L);%低通滤波后输出信号

        %调制
        [L1,g]=pulse(L,Fs,fd);         %将扩频信号转换为矩形波
        t=ts:ts:length(L);

        c1=sin(2*pi*Fc*t);                %载波信号
        y1=L1.*c1;                          %调制

        %加噪声
        y1=awgn(y1,snr);
        y1=y1+yg;

        %BPSK相干解调
        r1=y1;                  %接收到的BPSK信号
        cc1=sin(2*pi*Fc*t);     %本地载波
        rr=r1.*cc1;
        [B,A]=butter(4,wco);
        rrr=filtfilt(B,A,rr);   %低通滤波后输出信号

        %解扩
        mm=Fs/fd;
        %N=n/mm;
        rrrr=0;
        rrr0=zeros(112,1);
        rrr1=reshape(rrr,11200,N);
     
        for j=1:112
            rrr0(j,1)=rrr1(100*j-50,1);
        end
       
        rrr2=CM*rrr0(:,1);

        %抽样判决
        if rrr2>0
            rrrr=1;
         elseif rrr2<0
            rrrr=0;
        end

        %误码率计算
        if (rrrr-code2)~=0
            BER1=BER1+1;
        end
    end
    br(v,snr+1)=BER1;
end
v=v+1;
end
br1=br/NN;

