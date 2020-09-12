%直扩系统
clc;
clear;
f0=12;         %干扰频率
Fs=100;        %采样频率
Fc=10;        %载波频率
ts=1/Fs;       %采样间隔
fd=1;          %码元速率
N=100;         %码元个数
p=4;           %低通滤波器阶数

%低通滤波器截止频率
wc1=0.95;
wc2=0.9;
wc3=0.85;
wc4=0.6;
wc5=0.05;
gsr=10^(15/10);    %干信比
sgma=sqrt(2*N*gsr); %干扰能量

%产生N位信息码元
code2=randi([0,1],[1,N]);
code3=2.*code2-1;
[code1,m]=pulse(code3,Fs,fd);         %脉冲信号转换为矩形波
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
xlabel('HZ')
title('信源序列频谱');


%产生127位m序列
m=m_se([0 0 0 1 0 0 1]);         %利用移位寄存器产生m序列
m1=1-2.*m;
[m2,nn]=pulse(m1,Fs,fd);         %脉冲信号转换为矩形波
M=fft(m1)/length(m1);
M=shift(M,64,0);
subplot(323)
t=ts:ts:length(m1);
plot(t,m2);
axis([-inf,inf,-1.2,1.2])
title('m序列信号');
xlabel('n')
grid;
n=length(m1);
f=[0:Fs/n:Fs-Fs/n];
subplot(324)
plot(2*abs(fftshift(M)));
title('m序列频谱');

%生成扩频信号
L=kron(code3,m1);
subplot(325)
plot(L)
axis([0,127*10,-1.2,+1.2])
title('扩频后的序列');
xlabel('n')
[L1,g]=pulse(L,Fs,fd);         %将扩频信号转换为矩形波
M=fft(L1)/length(L1);
subplot(326)
n=length(L1);
f=[0:127/n:127-127/n]-127/2;
plot(f,2*abs(fftshift(M)));
title('扩频后序列频谱')
xlabel('HZ')
axis([0,+inf,-inf,+inf]);

[L1,m]=pulse(L,Fs,fd);         %将扩频信号转换为矩形波
t=ts:ts:length(L);

%画扩频后信号的特性图
figure(2)
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
L2=filtfilt(B,A,L1);%低通滤波后输出信号
figure(2)
subplot(423);
plot(t,L2);
axis([0,30,-1.2,1.2])
xlabel('n')
title('滤波平滑后的扩频信号（Wn=0.85）');
subplot(424)
M=fft(L2)/length(L2);
plot(f,2*abs(fftshift(M)));
title('滤波平滑后扩频信号频谱（Wn=0.85）')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

[B,A]=butter(p,wc4);
L2=filtfilt(B,A,L1);%低通滤波后输出信号
subplot(425);
plot(t,L2);
axis([0,30,-1.2,1.2])
xlabel('n')
title('滤波平滑后的扩频信号（Wn=0.8）');
subplot(426)
M=fft(L2)/length(L2);
plot(f,2*abs(fftshift(M)));
title('滤波平滑后扩频信号频谱（Wn=0.8）')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

[B,A]=butter(p,0.025);
L2=filtfilt(B,A,L1);%低通滤波后输出信号
subplot(427);
plot(t,L2);
axis([0,30,-1.2,1.2])
xlabel('n')
title('滤波平滑后的扩频信号（Wn=0.05）');
subplot(428)
M=fft(L2)/length(L2);
plot(f,2*abs(fftshift(M)));
title('滤波平滑后扩频信号频谱（Wn=0.05）')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')


%调制
c1=sin(2*pi*Fc*t);                %载波信号
y1=L2.*c1;                         

%产生单频干扰
yg=sgma*sin(2*pi*f0*t);

%画频谱
X1=fft(L2)/length(L2);
Y1=fft(y1)/length(y1);
C1=fft(c1)/length(c1);


figure(3)
subplot(321);
plot(t,L1);
axis([-inf,+inf,-1.2,1.2])
xlabel('n')
title('未调制的扩频信号');
grid;
subplot(322);
plot(f,2*abs(fftshift(X1)));
axis([0,+inf,-inf,+inf]);
xlabel('HZ')
title('未调制的扩频信号频谱');
grid;
subplot(323);
plot(c1);
axis([-inf,+inf,-1.2,1.2])
title('载波');
subplot(324);
plot(f+12,2*abs(fftshift(C1)));
title('载波频谱');
xlabel('HZ');


subplot(325);
plot(t,y1,t,L1,'-r');
axis([0,5,-1.2,1.2])
title('扩频后的BPSK信号');
grid;
subplot(326);
plot(f,2*abs(fftshift(Y1)));
axis([0,+inf,-inf,+inf]);
title('扩频后的BPSK信号频谱');
xlabel('HZ')
grid;

%加噪声
y2=awgn(y1,7);

%加单频干扰
%y1=y2+yg;


%BPSK相干解调
r1=y1;                   %接收到的BPSK信号
cc1=sin(2*pi*Fc*t);      %本地载波
rr=r1.*cc1;

%低通滤波
[B,A]=butter(2,0.1);
rrr=2*filtfilt(B,A,rr);   

%解扩
R=fft(r1)/length(r1);
R1=fft(y2)/length(y2);
RR=fft(rr)/length(rr);


%画解调后的信号特性
figure(8)
subplot(411)
plot(r1)
title('经过信道接收到的滤波平滑BPSK信号(前1000个点)')
xlabel('n')
axis([0,1000,-1.2,1.2])

%subplot(412)
%plot(f,abs(fftshift(R)));
%title('经过信道接收到的BPSK信号频谱')
%grid;

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

figure(5)
subplot(211)
plot(f,2*abs(fftshift(R1)));
title('经过信道接收到的BPSK信号频谱')
subplot(212)
plot(f,2*abs(fftshift(R)));
title('加入单频干扰的BPSK信号频谱')

%画解调后的信号特性
figure(5)
subplot(321)
plot(r1)
title('经过信道接收到的BPSK信号(snr=2)')
xlabel('n')
%axis([0,1000,-1.2,1.2])

subplot(322)
plot(f,2*abs(fftshift(R)));
axis([0,+inf,-inf,+inf])
title('经过信道接收到的BPSK信号频谱(snr=2)')
xlabel('HZ')
grid;

subplot(323)
plot(rr)
title('相干解调后信号')
xlabel('n')
%axis([0,1000,-1.2,1.2])

subplot(324)
plot(f,2*abs(fftshift(RR)));
axis([0,+inf,-inf,+inf])
title('相干解调后信号频谱')
xlabel('HZ')
grid;


mm=Fs/fd;
%N=n/mm;
rrrr=0;
rrr1=reshape(rrr,12700,N);
rrr2=zeros(1,N);
rrr0=zeros(127,N);

for i=1:N
    for j=1:127
         rrr0(j,i)=rrr1(100*j-50,i);
    end
end

%解扩
for i=1:N
    rrr2(i)=m1*rrr0(:,i);
end

subplot(325)
plot(rrr2)
title('解扩后信号')
xlabel('n')
axis([-inf,+inf,-1.2,1.2]);

[rrr3,g]=pulse(rrr2,Fs,fd);         %脉冲信号转换为矩形波
n=length(rrr3);
f=[0:1/n:1-1/n]-1/2;

RRR=fft(rrr3)/length(rrr3);
subplot(326)
plot(f,2*abs(fftshift(RRR)));
title('解扩后信号频谱')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')
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
