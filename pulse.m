function[x,m]=pulse(code,fs,fd);
m=fs/fd;
xx=[code;code];
xx=xx(:)';
x=xx;
for i=1:fix(m/2)-1;
    x=[x;xx];
end;
x=x(:)';