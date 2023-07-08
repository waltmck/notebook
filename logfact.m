function out=logfact(n)
persistent table;
if isempty(table), table=[0 cumsum(log(1:1000))]; end
neg=(n<0); n=max(0,n);
out=n;out(n<=1000)=table(n(n<=1000)+1);
x=n(n>1000); est=x.*log(x)-x+0.5*log(2*pi*x)+1./(12*x)-1./(360*x.^3);
out(n>1000)=est; out(neg)=inf;