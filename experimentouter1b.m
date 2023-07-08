function [val,grad]=experimentouter1b(u,al,be,yjs,hs,eps)
%attempt to depoissonize experimentouter1
if u<0 || u>1 || al<-1 || al>0, val=inf; grad=[nan nan nan]; return; end
i=(0:100)';
tp=1/u-1;
tm=sum(hs.*yjs);
fguess=@(k,tp,al)[.5-al*(2+tp)/4 k-(-k^2*(1+tp)*al).^(1/3) k+(-k^2*(1+tp)*al).^(1/3)];
for j=1:length(yjs)
    [x2,q,err]=inneredgesolveb(yjs(j),tp,al,be);
    if err==0, in(:,j)=[q;0;x2]; errs(:,j)=0; else [in(:,j),errs(:,j)]=bfgs2(@(in)experiment7c2(in,yjs(j),i,tp,al,be),fguess(yjs(j),tp,al)'); end%experimentinner1(yjs(j),u,al);
    [~,M,ders]=experiment7b2(in(1,j),in(2,j),in(3,j),yjs(j),i,tp,al,be);
    if in(2,j)==0, directs([1 3],:)=M([2 3],[1 3])\ders([2 3],:); else directs(:,:,j)=M\ders; end
end
if sum(errs)>1e-27, val=inf; grad=[nan nan nan]; return; end
qs=in(1,:); x1s=in(2,:); x2s=in(3,:);
%fobj=@(yjs,u,al,i,qs,x1s,x2s)(1-u).*(log(sum(exp(-u.*log((qs.*exp(logPois(x1s,i))+(1-qs).*exp(logPois(x2s,i)))./exp(logPois(yjs,i)))+logPois(x1s,i)),11))-al.*abs(x1s-yjs))+log(sum(exp((1-u).*log((qs.*exp(logPois(x1s,i))+(1-qs).*exp(logPois(x2s,i)))./exp(logPois(yjs,i)))+logPois(yjs,i)),11)).*u;
cis=log((qs.*exp(logPois(x1s,i))+(1-qs).*exp(logPois(x2s,i)))./exp(logPois(yjs,i)));
val=(eps*al+be*tm)*(1-u)+sum(hs.*((1-u).*(log(sum(exp(-u.*cis+logPois(x1s,i))))-al.*abs(x1s-yjs)-be*x1s)+log(sum(exp((1-u).*cis+logPois(yjs,i)))).*u));

poissonDeriv=@(L,k)exp(logPois(L,max(0,k-1))).*(1-2*(k==0)-(L./max(k,1)).*(k>0));
cidq=((exp(logPois(x1s,i))-exp(logPois(x2s,i)))./(qs.*exp(logPois(x1s,i))+(1-qs).*exp(logPois(x2s,i)))); %q derivative of each term of ci
cidx1=(qs.*poissonDeriv(x1s,i))./(qs.*exp(logPois(x1s,i))+(1-qs).*exp(logPois(x2s,i)));
cidx2=((1-qs).*poissonDeriv(x2s,i))./(qs.*exp(logPois(x1s,i))+(1-qs).*exp(logPois(x2s,i)));

deru=sum(hs.*((u-1).*sum(cis.*exp(-u*cis+logPois(x1s,i)),1)./sum(exp(-u*cis+logPois(x1s,i)),1)-u.*sum(cis.*exp((1-u)*cis+logPois(yjs,i)),1)./sum(exp((1-u)*cis+logPois(yjs,i)),1)-(log(sum(exp(-u.*cis+logPois(x1s,i))))-al.*abs(x1s-yjs)-be*x1s)+log(sum(exp((1-u).*cis+logPois(yjs,i))))),2);
derq=hs.*(sum(cidq.*exp(-u*cis+logPois(x1s,i)),1)./sum(exp(-u*cis+logPois(x1s,i)),1).*-u.*(1-u)+sum(cidq.*exp((1-u)*cis+logPois(yjs,i)),1)./sum(exp((1-u)*cis+logPois(yjs,i)),1).*u.*(1-u));
derx1=hs.*(sum(cidx1.*exp(-u*cis+logPois(x1s,i)),1)./sum(exp(-u*cis+logPois(x1s,i)),1).*-u.*(1-u)+sum(cidx1.*exp((1-u)*cis+logPois(yjs,i)),1)./sum(exp((1-u)*cis+logPois(yjs,i)),1).*u.*(1-u));
derx2=hs.*(sum(cidx2.*exp(-u*cis+logPois(x1s,i)),1)./sum(exp(-u*cis+logPois(x1s,i)),1).*-u.*(1-u)+sum(cidx2.*exp((1-u)*cis+logPois(yjs,i)),1)./sum(exp((1-u)*cis+logPois(yjs,i)),1).*u.*(1-u));
ders=[derq;derx1;derx2]; directs(~isfinite(directs))=0;%???
grad=-ders(:)'*reshape(permute(directs.*[-1/u^2 1 1],[1 3 2]),[],3)+[deru-eps*al-be*tm eps*(1-u)-sum(hs.*(1-u).*abs(x1s-yjs)) tm*(1-u)-sum(hs.*(1-u).*x1s)];

