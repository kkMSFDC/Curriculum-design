%扩频系统
clc;
clear;
f0=12;         %干扰频率
Fs=100;        %采样频率
Fc=10;         %载波频率
ts=1/Fs;       %采样间隔
fd=1;          %码元速率
N=100;         %码元个数
wc1=0.95;
wc2=0.9;
wc3=0.85;
wc4=0.8;
wc5=0.05;
p=4;
snr=7;
gsr=10^(15/10);    %能量比值
sgma=sqrt(2*N*gsr);

%产生信息码元
code2=randi([0,1],[1,N]);
code3=2.*code2-1;
[code1,g]=pulse(code3,Fs,fd);         %脉冲信号转换为矩形波
t=ts:ts:length(code3);
figure(1)
subplot(321)
plot(t,code1);
axis([-inf,+inf,-1.2,1.2])
title('信源信号');
xlabel('n')
grid;

n=length(code1);
f=[0:1/n:1-1/n]-1/2;
M=fft(code1)/length(code1);
subplot(322)
plot(f,2*abs(fftshift(M)));
axis([0,+inf,-inf,+inf]);
title('信源序列频谱');
xlabel('HZ')


%产生7位m序列
m=m_se([0 1 1]);
m1=1-2.*m;

%产生16位Walsh序列

wals=walsh_se(16);
wals1=wals(8,:);

[m2,nn]=pulse(wals1,Fs,fd);         %脉冲信号转换为矩形波
M=fft(wals1)/length(wals1);
M=shift(M,-8,0);
subplot(323)
t=ts:ts:length(wals1);
plot(t,m2);
axis([0,inf,-1.2,1.2])
xlabel('n')
title('walsh序列信号');
grid;
n=length(wals1);
f=(0:16/n:16-16/n)-16/2;
subplot(324)
plot(f,2*abs(fftshift(M)));
axis([0,+inf,-inf,+inf])
title('walsh序列频谱');
figure(6)
Rg=xg(wals1,wals1);
plot(0:16,Rg)
title('walsh序列的自相关函数')
xlabel('偏移量')

%产生复合序列
CM=kron(wals1,m1);

%扩频
L=kron(code3,CM);
figure(1)
subplot(325)
stem(L)
axis([0,127*10,-1.2,+1.2])
title('扩频序列信号');
[L1,g]=pulse(L,Fs,fd);         %将扩频信号转换为矩形波
M=fft(L1)/length(L1);
subplot(326)
n=length(L1);
f=[0:112/n:112-112/n]-112/2;
plot(f,2*abs(fftshift(M)));
title('扩频序列频谱')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')



%调制
[L1,m]=pulse(L,Fs,fd);         %将扩频信号转换为矩形波
t=ts:ts:length(L);

c1=sin(2*pi*Fc*t);                %载波信号
y1=L1.*c1;                          %调制


figure(9)
subplot(421);
plot(t,L1);
axis([0,30,-1.2,1.2])
xlabel('n')
title('未滤波平滑的扩频信号');
subplot(422)
M=fft(L1)/length(L1);
plot(f,2*abs(fftshift(M)));
title('未滤波平滑的扩频信号频谱')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

%平滑
[B,A]=butter(p,wc3);
L1=filtfilt(B,A,L1);%低通滤波后输出信号
figure(9)
subplot(423);
plot(t,L1);
axis([0,30,-1.2,1.2])
xlabel('n')
title('滤波平滑后的扩频信号（Wn=0.85）');
subplot(424)
M=fft(L1)/length(L1);
plot(f,2*abs(fftshift(M)));
title('滤波平滑后扩频信号频谱（Wn=0.85）')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

[B,A]=butter(p,wc4);
L1=filtfilt(B,A,L1);%低通滤波后输出信号
subplot(425);
plot(t,L1);
axis([0,30,-1.2,1.2])
xlabel('n')
title('滤波平滑后的扩频信号（Wn=0.8）');
subplot(426)
M=fft(L1)/length(L1);
plot(f,2*abs(fftshift(M)));
title('滤波平滑后扩频信号频谱（Wn=0.8）')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

[B,A]=butter(p,0.05);
L1=filtfilt(B,A,L1);%低通滤波后输出信号
subplot(427);
plot(t,L1);
axis([0,30,-1.2,1.2])
xlabel('n')
title('滤波平滑后的扩频信号（Wn=0.05）');
subplot(428)
M=fft(L1)/length(L1);
plot(f,2*abs(fftshift(M)));
title('滤波平滑后扩频信号频谱（Wn=0.05）')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

%产生单频干扰
yg=sgma*sin(2*pi*f0*t);

%画频谱
X1=fft(L1)/length(L1);
Y1=fft(y1)/length(y1);
C1=fft(c1)/length(c1);

n=length(L1);
f=[0:Fs/n:Fs-Fs/n]-Fs/2;

figure(2)
subplot(321);
plot(t,L1);
axis([-inf,+inf,-1.2,1.2])
xlabel('n')
title('未调制扩频信号');
grid;
subplot(322);
plot(f,2*abs(fftshift(X1)));
axis([0,+inf,-inf,+inf])
xlabel('HZ')
title('未调制的扩频信号频谱');
grid;

subplot(323);
plot(c1);
axis([-inf,+inf,-1.2,1.2])
title('载波');
xlabel('n')
grid;
subplot(324);
plot(f+10,2*abs(fftshift(C1)));
axis([-10,+inf,-inf,+inf])
xlabel('HZ')
title('载波频谱');
grid;

subplot(325);
plot(t,y1,t,L1,'-r');
axis([0,10,-1.2,1.2])
title('BPSK信号');
grid;
subplot(326);
plot(f,2*abs(fftshift(Y1)));
axis([0,+inf,-inf,+inf])
xlabel('HZ')
title('扩频后的BPSK信号频谱');
grid;

%加噪声
y1=awgn(y1,snr);

%加单频干扰
y1=y1+yg;

%BPSK相干解调
r1=y1;                  %接收到的BPSK信号
cc1=sin(2*pi*Fc*t);     %本地载波
rr=r1.*cc1;
[B,A]=butter(2,0.1);
rrr=2*filtfilt(B,A,rr);   %低通滤波后输出信号

%解扩
R=fft(r1)/length(r1);
R1=fft(y1)/length(y1);
RR=fft(rr)/length(rr);


%画解调后的信号特性
figure(8)
subplot(411)
plot(r1)
title('经过信道接收到的滤波平滑后BPSK信号(前1000个点)')
xlabel('n')
axis([0,1000,-1.2,1.2])

subplot(412)
plot(rr)
title('相干解调后信号')
axis([0,1000,-1.2,1.2])
xlabel('n')

subplot(413)
plot(rrr)
title('经过低通滤波后信号')
axis([0,1000,-1.2,1.2])
xlabel('n')

%画解调后的信号特性
figure(3)
subplot(321)
plot(r1)
title('经过信道接收到的BPSK信号(前1000个点)')
xlabel('n')
%axis([0,1000,-1.2,1.2])

subplot(322)
plot(f,2*abs(fftshift(R)));
axis([0,+inf,-inf,+inf])
xlabel('HZ')
title('经过信道接收到的BPSK信号频谱')
grid;

subplot(323)
plot(rr)
title('相干解调后信号')
xlabel('n')
axis([0,1000,-1.2,1.2])

subplot(324)
plot(f,2*abs(fftshift(RR)));
axis([0,+inf,-inf,+inf])
xlabel('HZ')
title('相干解调后信号频谱')
grid;


mm=Fs/fd;
rrrr=0;
rrr1=reshape(rrr,11200,N);
rrr2=zeros(1,N);
rrr0=zeros(112,N);

for i=1:N
    for j=1:112
         rrr0(j,i)=rrr1(100*j-50,i);
    end
end

%解扩
for i=1:N
    rrr2(i)=CM*rrr0(:,i);
end

subplot(325)
plot(rrr2)
title('解扩后信号')
xlabel('n')
axis([-inf,+inf,-1.2,1.2]);

[rrr3,g]=pulse(rrr2,Fs,fd);         %脉冲信号转换为矩形波
n=length(rrr3);
f=[0:Fs/n:Fs-Fs/n]-Fs/2;

RRR=fft(rrr3)/length(rrr3);
subplot(326)
plot(f,2*abs(fftshift(RRR)));
xlabel('HZ')
title('解扩后信号频谱')
axis([0,+inf,-inf,+inf]);
grid;


%抽样判决
for k=1:N
    if rrr2(k)>0
        rrrr(k)=1;
    elseif rrr2(k)<0
        rrrr(k)=0;
    end
end
figure(8)
subplot(414)
[yy,m]=pulse(rrrr,100,fd);         %脉冲信号转换为矩形波
N=length(rrrr);
m=linspace(0,N,N*100);
plot(m,yy)
title('输出信号')
xlabel('n')
axis([0,100,-0.2,1.2])

%误码率计算
BER=0;
for i=1:N
    if (rrrr(i)-code2(i))==0
        BER=BER+1;
    end
end
BER=N-BER;
BER=BER/length(code2)
