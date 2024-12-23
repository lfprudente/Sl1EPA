function [nl,ni] = evalnlni(n,m,x,lambda,lambdaesc,g,sc,scaling)

ni = zeros(n,1);
nl = g;

for i = 1:m
    [~,ncs,~] = sevalnc(n,x,i,sc,scaling);
    ni = ni + lambdaesc(i) * ncs;
    nl = nl + lambda(i) * ncs;
        
end

[nl] = reshapevector(nl);
[ni] = reshapevector(ni);