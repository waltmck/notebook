function [v,M,ders]=experiment7(q,x1,x2,yj,i,tp,al,be)
%like experiment7b but with another dual variable, beta, for
%depoissonization

%computes 3 errors, and 3x3 derivative matrix, for optimization/solving;
%see also experiment7c that makes this into a quadratic objective function.
%3rd output argument is derivatives with respect to tp and al

% poissonDeriv=@(L,k)exp(logPois(L,max(0,k-1))).*(1-2*(k==0)-(L./max(k,1)).*(k>0));
% poissonDerivxx=@(L,k)exp(logPois(L,max(0,k-2))).*(k.^2-k.*(2*L+1)+L.^2)./(max(1,k).*max(1,k-1)).*(k>1)+exp(-L).*((k==0)+(L-2).*(k==1));
%
% lpoix1=logPois(x1,i); poix1=exp(lpoix1);
% poidx1=poissonDeriv(x1,i);
% poiddx1=poissonDerivxx(x1,i);
% lpoix2=logPois(x2,i); poix2=exp(lpoix2);
% poidx2=poissonDeriv(x2,i);
% poiddx2=poissonDerivxx(x2,i);
% lpoik=logPois(k,i);

fi=1./cumprod([1;(1:length(i)-1)']);
im2=i(3:end);
poix1=exp(-x1+i.*log(x1+(x1==0 & i==0))).*fi;
poidx1=[-poix1(1); poix1(1:end-1).*(1-(x1./i(2:end)))];
poiddx1=[poix1(1)*[1; x1-2]; poix1(1:end-2).*(im2.^2-im2.*(2*x1+1)+x1.^2)./(im2.*(im2-1))];
poix2=exp(-x2+i*log(x2)).*fi;
poidx2=[-poix2(1); poix2(1:end-1).*(1-(x2./i(2:end)))];
poiddx2=[poix2(1)*[1; x2-2]; poix2(1:end-2).*(im2.^2-im2.*(2*x2+1)+x2.^2)./(im2.*(im2-1))];
poik=exp(-yj+i*log(yj)).*fi;

ci=log((poix1*q+poix2*(1-q))./poik); cip=exp(-ci/(1+tp));
cidx1=q.*poidx1./(poix1*q+poix2*(1-q));
cidx2=(1-q).*poidx2./(poix1*q+poix2*(1-q));
cidp=(poix1-poix2)./(poix1*q+poix2*(1-q));
poix12=[poix1 poix2]; temp=1./(sum(cip.*poix12)*(1+tp));
valdcix1=-sum(cidx1.*cip.*poix12).*temp;
valdcix2=-sum(cidx2.*cip.*poix12).*temp;
valdcip=-sum(cidp.*cip.*poix12).*temp;
derx1=sum(cip.*poidx1)./sum(cip.*poix1);
derx2=sum(cip.*poidx2)./sum(cip.*poix2);

derx1dcix1=(sum(cip.*poidx1).*sum(cidx1.*cip.*poix1)-sum(cidx1.*cip.*poidx1).*sum(cip.*poix1))./sum(cip.*poix1).^2/(1+tp);
derx1dcix2=(sum(cip.*poidx1).*sum(cidx2.*cip.*poix1)-sum(cidx2.*cip.*poidx1).*sum(cip.*poix1))./sum(cip.*poix1).^2/(1+tp);
derx1dcip=(sum(cip.*poidx1).*sum(cidp.*cip.*poix1)-sum(cidp.*cip.*poidx1).*sum(cip.*poix1))./sum(cip.*poix1).^2/(1+tp);
derx1ddx1=(sum(cip.*poix1).*sum(cip.*poiddx1)-sum(cip.*poidx1).^2)./sum(cip.*poix1).^2;

derx2dcix1=(sum(cip.*poidx2).*sum(cidx1.*cip.*poix2)-sum(cidx1.*cip.*poidx2).*sum(cip.*poix2))./sum(cip.*poix2).^2/(1+tp);
derx2dcix2=(sum(cip.*poidx2).*sum(cidx2.*cip.*poix2)-sum(cidx2.*cip.*poidx2).*sum(cip.*poix2))./sum(cip.*poix2).^2/(1+tp);
derx2dcip=(sum(cip.*poidx2).*sum(cidp.*cip.*poix2)-sum(cidp.*cip.*poidx2).*sum(cip.*poix2))./sum(cip.*poix2).^2/(1+tp);
derx2ddx2=(sum(cip.*poix2).*sum(cip.*poiddx2)-sum(cip.*poidx2).^2)./sum(cip.*poix2).^2;

err1dx1=derx1dcix1+derx1ddx1;
err1dx2=derx1dcix2;
err1dp=derx1dcip;
err2dx1=derx2dcix1;
err2dx2=derx2dcix2+derx2ddx2;
err2dp=derx2dcip;
err3dx1=diff(valdcix1)-derx1+al*sign(x1-yj)+be;
err3dx2=diff(valdcix2)+derx2-al*sign(x2-yj)-be;
err3dp=diff(valdcip);
M=[err1dp err1dx1 err1dx2; err2dp err2dx1 err2dx2; err3dp err3dx1 err3dx2];

err1=sum(cip.*poidx1)./sum(cip.*poix1)+al-be;
err2=sum(cip.*poidx2)./sum(cip.*poix2)-al-be;
err3=diff(log(sum(cip.*poix12))-al*abs([x1 x2]-yj)-be*[x1 x2]);
v=[err1 err2 err3]; %if q<=0 || x1>k || x2<k, v(:)=1+2i;end

if nargout>=3,
    cipdtp=cip.*ci/(1+tp)^2;
    err1dtp=sum(cipdtp.*poidx1)./sum(cip.*poix1)-sum(cip.*poidx1)*sum(cipdtp.*poix1)/sum(cip.*poix1)^2;
    err2dtp=sum(cipdtp.*poidx2)./sum(cip.*poix2)-sum(cip.*poidx2)*sum(cipdtp.*poix2)/sum(cip.*poix2)^2;
    err3dtp=diff(sum(cipdtp.*poix12)./sum(cip.*poix12));
    ders=[[err1dtp;err2dtp;err3dtp] [1;-1;-(x1+x2-2*yj)] [-1;-1;x1-x2]];
end