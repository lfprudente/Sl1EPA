function [nal,store] = evalnal(x,store,n,m,equatn,tau,rho,sc,scaling)

global evalfun evalgrad nevalnal typephi

nevalnal = nevalnal + 1;

[~,gs,flag] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

if ~isfield(store, 'c')
    cs = [];
    for i = 1:m
        [~,cs(i),flag] = sevalc(n,x,i,sc,scaling);
        evalfun = evalfun + 1;
    end
    cs = cs';

    store.c = cs;       
end

cs = store.c;

nal = [];

nal = gs;
for i = 1:m
    [dphival] = dphi(cs(i),tau,typephi);

    [~,ncs,flag] = sevalnc(n,x,i,sc,scaling);
    evalgrad = evalgrad + 1;

    if ( equatn(i) == true )
        lambda = rho(i) * dphival;
    else
        lambda = rho(i) * 0.5 * ( 1 + dphival );
        % lambda = rho(i) * ( 1 + dphival );
    end
    
    nal = nal + lambda * ncs;
end

[nal] = reshapevector(nal);