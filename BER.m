%��Ƶϵͳ
br=zeros(2,6);
NN=100000;
v=1;
wco=0.5;
%����127λm����
m=m_se([0 0 0 1 0 0 1]);
m1=1-2*m;

%snr=0;
%for wc=0.1:0.2:0.3;
for snr=0:2:10
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
         L=kron(code3,m1);

        %����
        [L1,g]=pulse(L,Fs,fd);         %����Ƶ�ź�ת��Ϊ���β�
        t=ts:ts:length(L);

        c1=sin(2*pi*Fc*t);                %�ز��ź�
        y1=L1.*c1;                          %����

        %������
        y1=awgn(y1,snr);

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
        rrr0=zeros(127,1);
        rrr1=reshape(rrr,12700,N);
     
        for j=1:127
            rrr0(j,1)=rrr1(100*j-50,1);
        end
       
        rrr2=m1*rrr0(:,1);

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
    br(v,snr/2+1)=BER1;
end
%v=v+1;
%end
br1=br/NN;

