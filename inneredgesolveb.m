function [x2,q,err]=inneredgesolve(yj,tp,al,be)
tt=-1/(1+tp);
x2=1;
for i=1:50
    ff=expm1(-x2*(x2/yj)^tt)/expm1(-(x2+al*(x2-2*yj)+be*x2))*(yj/x2)^tt-1/(1+al+be);
    if ff>0, if i==1, L=0; else L=x2*.5; end; break; end
    x2=x2*2;
end
if i==50, x2=nan; q=nan; err=1; return; end
U=x2;
for i=1:50
    x2=(L+U)/2;
    ff=expm1(-x2*(x2/yj)^tt)/expm1(-(x2+al*(x2-2*yj)+be*x2))*(yj/x2)^tt-1/(1+al+be);
    if ff>0, U=x2; else L=x2; end
end
q=1./(1+1./(exp(-x2).*(-1+((exp(x2.*(x2./yj).^tt)-1)./(exp(x2+al.*(x2-2*yj)+be*x2)-1)).^(1/tt))));
test=x2+al*(x2-2*yj)+be*x2-x2*(x2/yj)^tt;
if test<0 || test>log((1-al+be)/(1+al+be)), err=1; disp(['test error: ' num2str(test)]); else err=0; end

