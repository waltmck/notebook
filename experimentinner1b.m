function [in,err]=experimentinner1(yj,u,al,be,in)
%one issue with this code is when q converges to <1e-15, then pinv(M) is
%doing the wrong thing
tp=1/u-1;
if al==0, in=[.5 yj yj]; err=0; return; end %check
if nargin<5, [x2,q,err]=inneredgesolveb(yj,tp,al,be); in=[q yj*.1 x2];  if err==0, in(2)=0; return; end; end
if any(~isfinite(in)) || any(abs(in)<1e-100),
    if -al*(2+tp)/4<.25 && (-yj^2*(1+tp)*al)^(1/3)<.75*yj,
        in=[.5-al*(2+tp)/4 yj-(-yj^2*(1+tp)*al).^(1/3) yj+(-yj^2*(1+tp)*al).^(1/3)];
    else
        in=[.5 .5*yj 1.5*yj];
    end;
end
i=(0:100)';
lb=[0 0 yj]; ub=[1 yj inf]; th=.5;
for j=1:50,
    [v,M]=experiment7b2(in(1),in(2),in(3),yj,i,tp,al,be); if sum(v.^2)<1e-30, break; end
    p=-v/M';%-(pinv(M)*v')';%-v/M';%(inv(M)*v')'; %I forget when M will be singular
    pu=(ub-in)*th; pl=(lb-in)*th;
    a=min(1,1./max(max(p./[pl;pu])));
    if a<.01,
        [m,k]=max(p./pl); if m>=100, xt=in(k)*(1-th)+lb(k)*th; ind=k; end; [m,k]=max(p./pu); if m>=100, xt=in(k)*(1-th)+ub(k)*th; ind=k; end;
        temp=pinv(M(:,[1:ind-1 ind+1:end]))*(v+M(:,ind)'*(xt-in(ind)))';
        p=zeros(1,3);p([1:ind-1 ind+1:end])=-temp;p(ind)=xt-in(ind);
        a=min(1,1./max(max(p./[pl;pu])));
        if a<.01, disp([j,a]), end
    end
    in=in+p*a;
end
[v,M]=experiment7b2(in(1),in(2),in(3),yj,i,tp,al,be); err=sum(v.^2);