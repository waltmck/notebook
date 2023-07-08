function M=logPois(lambda,k) 
% real vector x=x1,x2,x3,...,xk, and integer vector n=n1, n2,...,nm
% output matrix M of poissons, where M(i,j)=pois(xi,nj)=Probability of a
% poisson with expectation xi is equal to nj.
%in a typical run, maybe call makePois(x*k,1:r)

%big entries do the stirling approx of factorial
%[X,Y]=meshgrid(x,n);
%M=(log(X).*Y-X-logfact(Y))';
M=(log(lambda).*k-lambda-logfact(k));%singleton expansion - newer matlabs
M((lambda==0) & (k==0))=0;