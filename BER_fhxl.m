%��Ƶϵͳ
br=zeros(9,11);
NN=100000;
N=1;
v=1;
wco=0.5;
wc=0.8;
p=4;
Fs=100;
ts=1/Fs;

%����7λm����
m=m_se([0 1 1]);
m1=1-2*m;

%����16λWalsh����
wals=walsh_se(16);
wals1=wals(8,:);

%����112λ��������
CM=kron(wals1,m1);

gsr=10^(15/10);    %������ֵ
sgma=sqrt(2*N*gsr);
t=ts:ts:112;
%������Ƶ����

%snr=0;
for f0=10:1:10
    yg=sgma*sin(2*pi*f0*t);
for snr=0:1:10
    BER1=0;
    for i=1:NN
        Fs=100;        %����Ƶ��
        Fc=10;         %�ز�Ƶ��
        ts=1/Fs;       %�������
        fd=1;          %��Ԫ����
        N=1;         %��Ԫ����

        %������Ϣ��Ԫ
        code2=randi([0,1],[1,N]);
        code3=2.*code2-1;

       
        %[m2,nn]=pulse(m1,Fs,fd);         %�����ź�ת��Ϊ���β�

         %��Ƶ
         L=kron(code3,CM);
         
         %ƽ��
         [B,A]=butter(p,wc);
         L=filtfilt(B,A,L);%��ͨ�˲�������ź�

        %����
        [L1,g]=pulse(L,Fs,fd);         %����Ƶ�ź�ת��Ϊ���β�
        t=ts:ts:length(L);

        c1=sin(2*pi*Fc*t);                %�ز��ź�
        y1=L1.*c1;                          %����

        %������
        y1=awgn(y1,snr);
        y1=y1+yg;

        %BPSK��ɽ��
        r1=y1;                  %���յ���BPSK�ź�
        cc1=sin(2*pi*Fc*t);     %�����ز�
        rr=r1.*cc1;
        [B,A]=butter(4,wco);
        rrr=filtfilt(B,A,rr);   %��ͨ�˲�������ź�

        %����
        mm=Fs/fd;
        %N=n/mm;
        rrrr=0;
        rrr0=zeros(112,1);
        rrr1=reshape(rrr,11200,N);
     
        for j=1:112
            rrr0(j,1)=rrr1(100*j-50,1);
        end
       
        rrr2=CM*rrr0(:,1);

        %�����о�
        if rrr2>0
            rrrr=1;
         elseif rrr2<0
            rrrr=0;
        end

        %�����ʼ���
        if (rrrr-code2)~=0
            BER1=BER1+1;
        end
    end
    br(v,snr+1)=BER1;
end
v=v+1;
end
br1=br/NN;

