function [fval,grad]=experiment7(in,k,i,tp,al,be)
[v,M]=experiment7b2(in(1),in(2),in(3),k,i,tp,al,be);
fval=sum(v.^2);
grad=M'*v'*2;