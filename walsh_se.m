function walsh = walsh_se(N)
M=ceil(log2(N));
Hn=1;
for n=1:M
    H2n=[Hn,Hn;Hn,-Hn];
    Hn=H2n;
end
walsh=Hn;
end
