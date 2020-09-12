%ֱ��ϵͳ
clc;
clear;
f0=12;         %����Ƶ��
Fs=100;        %����Ƶ��
Fc=10;        %�ز�Ƶ��
ts=1/Fs;       %�������
fd=1;          %��Ԫ����
N=100;         %��Ԫ����
p=4;           %��ͨ�˲�������

%��ͨ�˲�����ֹƵ��
wc1=0.95;
wc2=0.9;
wc3=0.85;
wc4=0.6;
wc5=0.05;
gsr=10^(15/10);    %���ű�
sgma=sqrt(2*N*gsr); %��������

%����Nλ��Ϣ��Ԫ
code2=randi([0,1],[1,N]);
code3=2.*code2-1;
[code1,m]=pulse(code3,Fs,fd);         %�����ź�ת��Ϊ���β�
t=ts:ts:length(code3);
figure(1)
subplot(321)
plot(t,code1);
axis([-inf,+inf,-1.2,1.2])
title('��Դ�ź�');
xlabel('n')
grid;
n=length(code1);
f=[0:1/n:1-1/n]-1/2;
M=fft(code1)/length(code1);
subplot(322)
plot(f,2*abs(fftshift(M)));
axis([0,+inf,-inf,+inf]);
xlabel('HZ')
title('��Դ����Ƶ��');


%����127λm����
m=m_se([0 0 0 1 0 0 1]);         %������λ�Ĵ�������m����
m1=1-2.*m;
[m2,nn]=pulse(m1,Fs,fd);         %�����ź�ת��Ϊ���β�
M=fft(m1)/length(m1);
M=shift(M,64,0);
subplot(323)
t=ts:ts:length(m1);
plot(t,m2);
axis([-inf,inf,-1.2,1.2])
title('m�����ź�');
xlabel('n')
grid;
n=length(m1);
f=[0:Fs/n:Fs-Fs/n];
subplot(324)
plot(2*abs(fftshift(M)));
title('m����Ƶ��');

%������Ƶ�ź�
L=kron(code3,m1);
subplot(325)
plot(L)
axis([0,127*10,-1.2,+1.2])
title('��Ƶ�������');
xlabel('n')
[L1,g]=pulse(L,Fs,fd);         %����Ƶ�ź�ת��Ϊ���β�
M=fft(L1)/length(L1);
subplot(326)
n=length(L1);
f=[0:127/n:127-127/n]-127/2;
plot(f,2*abs(fftshift(M)));
title('��Ƶ������Ƶ��')
xlabel('HZ')
axis([0,+inf,-inf,+inf]);

[L1,m]=pulse(L,Fs,fd);         %����Ƶ�ź�ת��Ϊ���β�
t=ts:ts:length(L);

%����Ƶ���źŵ�����ͼ
figure(2)
subplot(421);
plot(t,L1);
axis([0,30,-1.2,1.2])
xlabel('n')
title('δ�˲�ƽ������Ƶ�ź�');
subplot(422)
M=fft(L1)/length(L1);
plot(f,2*abs(fftshift(M)));
title('δ�˲�ƽ������Ƶ�ź�Ƶ��')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

%ƽ��
[B,A]=butter(p,wc3);
L2=filtfilt(B,A,L1);%��ͨ�˲�������ź�
figure(2)
subplot(423);
plot(t,L2);
axis([0,30,-1.2,1.2])
xlabel('n')
title('�˲�ƽ�������Ƶ�źţ�Wn=0.85��');
subplot(424)
M=fft(L2)/length(L2);
plot(f,2*abs(fftshift(M)));
title('�˲�ƽ������Ƶ�ź�Ƶ�ף�Wn=0.85��')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

[B,A]=butter(p,wc4);
L2=filtfilt(B,A,L1);%��ͨ�˲�������ź�
subplot(425);
plot(t,L2);
axis([0,30,-1.2,1.2])
xlabel('n')
title('�˲�ƽ�������Ƶ�źţ�Wn=0.8��');
subplot(426)
M=fft(L2)/length(L2);
plot(f,2*abs(fftshift(M)));
title('�˲�ƽ������Ƶ�ź�Ƶ�ף�Wn=0.8��')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')

[B,A]=butter(p,0.025);
L2=filtfilt(B,A,L1);%��ͨ�˲�������ź�
subplot(427);
plot(t,L2);
axis([0,30,-1.2,1.2])
xlabel('n')
title('�˲�ƽ�������Ƶ�źţ�Wn=0.05��');
subplot(428)
M=fft(L2)/length(L2);
plot(f,2*abs(fftshift(M)));
title('�˲�ƽ������Ƶ�ź�Ƶ�ף�Wn=0.05��')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')


%����
c1=sin(2*pi*Fc*t);                %�ز��ź�
y1=L2.*c1;                         

%������Ƶ����
yg=sgma*sin(2*pi*f0*t);

%��Ƶ��
X1=fft(L2)/length(L2);
Y1=fft(y1)/length(y1);
C1=fft(c1)/length(c1);


figure(3)
subplot(321);
plot(t,L1);
axis([-inf,+inf,-1.2,1.2])
xlabel('n')
title('δ���Ƶ���Ƶ�ź�');
grid;
subplot(322);
plot(f,2*abs(fftshift(X1)));
axis([0,+inf,-inf,+inf]);
xlabel('HZ')
title('δ���Ƶ���Ƶ�ź�Ƶ��');
grid;
subplot(323);
plot(c1);
axis([-inf,+inf,-1.2,1.2])
title('�ز�');
subplot(324);
plot(f+12,2*abs(fftshift(C1)));
title('�ز�Ƶ��');
xlabel('HZ');


subplot(325);
plot(t,y1,t,L1,'-r');
axis([0,5,-1.2,1.2])
title('��Ƶ���BPSK�ź�');
grid;
subplot(326);
plot(f,2*abs(fftshift(Y1)));
axis([0,+inf,-inf,+inf]);
title('��Ƶ���BPSK�ź�Ƶ��');
xlabel('HZ')
grid;

%������
y2=awgn(y1,7);

%�ӵ�Ƶ����
%y1=y2+yg;


%BPSK��ɽ��
r1=y1;                   %���յ���BPSK�ź�
cc1=sin(2*pi*Fc*t);      %�����ز�
rr=r1.*cc1;

%��ͨ�˲�
[B,A]=butter(2,0.1);
rrr=2*filtfilt(B,A,rr);   

%����
R=fft(r1)/length(r1);
R1=fft(y2)/length(y2);
RR=fft(rr)/length(rr);


%���������ź�����
figure(8)
subplot(411)
plot(r1)
title('�����ŵ����յ����˲�ƽ��BPSK�ź�(ǰ1000����)')
xlabel('n')
axis([0,1000,-1.2,1.2])

%subplot(412)
%plot(f,abs(fftshift(R)));
%title('�����ŵ����յ���BPSK�ź�Ƶ��')
%grid;

subplot(412)
plot(rr)
title('��ɽ�����ź�')
axis([0,1000,-1.2,1.2])
xlabel('n')

subplot(413)
plot(rrr)
title('������ͨ�˲����ź�')
axis([0,1000,-1.2,1.2])
xlabel('n')

figure(5)
subplot(211)
plot(f,2*abs(fftshift(R1)));
title('�����ŵ����յ���BPSK�ź�Ƶ��')
subplot(212)
plot(f,2*abs(fftshift(R)));
title('���뵥Ƶ���ŵ�BPSK�ź�Ƶ��')

%���������ź�����
figure(5)
subplot(321)
plot(r1)
title('�����ŵ����յ���BPSK�ź�(snr=2)')
xlabel('n')
%axis([0,1000,-1.2,1.2])

subplot(322)
plot(f,2*abs(fftshift(R)));
axis([0,+inf,-inf,+inf])
title('�����ŵ����յ���BPSK�ź�Ƶ��(snr=2)')
xlabel('HZ')
grid;

subplot(323)
plot(rr)
title('��ɽ�����ź�')
xlabel('n')
%axis([0,1000,-1.2,1.2])

subplot(324)
plot(f,2*abs(fftshift(RR)));
axis([0,+inf,-inf,+inf])
title('��ɽ�����ź�Ƶ��')
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

%����
for i=1:N
    rrr2(i)=m1*rrr0(:,i);
end

subplot(325)
plot(rrr2)
title('�������ź�')
xlabel('n')
axis([-inf,+inf,-1.2,1.2]);

[rrr3,g]=pulse(rrr2,Fs,fd);         %�����ź�ת��Ϊ���β�
n=length(rrr3);
f=[0:1/n:1-1/n]-1/2;

RRR=fft(rrr3)/length(rrr3);
subplot(326)
plot(f,2*abs(fftshift(RRR)));
title('�������ź�Ƶ��')
axis([0,+inf,-inf,+inf]);
xlabel('HZ')
grid;


%�����о�
for k=1:N
    if rrr2(k)>0
        rrrr(k)=1;
    elseif rrr2(k)<0
        rrrr(k)=0;
    end
end

figure(8)
subplot(414)
[yy,m]=pulse(rrrr,100,fd);         %�����ź�ת��Ϊ���β�
N=length(rrrr);
m=linspace(0,N,N*100);
plot(m,yy)
title('����ź�')
xlabel('n')
axis([0,100,-0.2,1.2])

%�����ʼ���
BER=0;
for i=1:N
    if (rrrr(i)-code2(i))==0
        BER=BER+1;
    end
end
BER=N-BER;
BER=BER/length(code2)
