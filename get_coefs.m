function cis=get_coefs(n, k, eps)
    i=(0:k)';
    yjs=[k/n]; hs=[1]; eps=eps*k/n;
    [in,fval]=bfgs2(@(in)experimentouter1b(in(1),in(2),in(3),yjs,hs,eps),[.5 0 0]');
    [val,grad]=experimentouter1b(in(1),in(2),in(3),yjs,hs,eps);
    fguess=@(k,tp,al)[.5-al*(2+tp)/4 k-(-k^2*(1+tp)*al).^(1/3) k+(-k^2*(1+tp)*al).^(1/3)];
    [qx1x2,err]=bfgs2(@(in2)experiment7c2(in2,yjs,i,1/in(1)-1,in(2),in(3)),fguess(yjs,1/in(1)-1,in(2))');
    cis=log((qx1x2(1).*exp(logPois(qx1x2(2),i))+(1-qx1x2(1)).*exp(logPois(qx1x2(3),i)))./exp(logPois(yjs,i)));
end
