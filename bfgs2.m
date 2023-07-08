function [x,fval]=bfgs2(f,x)
%my attempt at bfgs minimization, *assuming* the global min has function
%value 0 - which is weird in general; here it's designed for objective
%function that is sum of squares of constraint violation.
H=eye(length(x));
[fval,grad]=f(x);
c1=1e-4; c2=.9; fe=0;%total function evaluations, out of curiosity
reset=0; resets=0; maxresets=1;%5;
%loop:
for j=1:450
    p=-H*grad(:);
    %line search here!
    a=1; i=0;
    [newfval,newgrad]=f(x+a*p); %if sum(newgrad.^2)<1e-29, return; end%???
    if isreal(newfval) && newfval<=fval+c1*a*p'*grad(:) && p'*newgrad(:)<1e-28, mul=2; else mul=.5; end
    if mul==.5 || abs(p'*newgrad(:))>c2*abs(p'*grad(:)),
        while 1
            i=i+1; a=a*mul; [newfval,newgrad]=f(x+a*p);
            if isreal(newfval) && newfval<=fval+c1*a*p'*grad(:) && abs(p'*newgrad(:))<=c2*abs(p'*grad(:)), break; end
            if mul>1 && p'*newgrad(:)>0, a=a/mul; [newfval,newgrad]=f(x+a*p); disp('special case!'); break; end
            if i==30, 
                H=eye(length(x));reset=1; break;%return; 
            end %can't find a point that satisfies the Wolfe conditions
        end
    end
    if any(isnan(newgrad)), return; end%???
    fe=fe+1+i;
    s=a*p;
    x=x+s;
    oldgrad=grad;
    fval=newfval; grad=newgrad; %[fval,grad]=f(x);
%    if sum(grad.^2)<1e-30, return; end%???
    y=grad(:)-oldgrad(:); if sum(y.^2)<1e-30, return; end;%???
    if reset, reset=0; resets=resets+1; if resets==maxresets, return; end; else H=H+(s'*y+y'*H*y)*(s*s')/(s'*y)^2-(H*y*s'+s*y'*H)/(s'*y); end
end