function [al,store] = evalal(x,store,n,m,equatn,tau,rho,sc,scaling)

global evalfun nevalal typephi

nevalal = nevalal + 1;

[~,fs,flag] = sevalf(n,x,sc,scaling);
evalfun = evalfun + 1;

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

al = fs;
for i = 1:m
    [phival] = phi(cs(i),tau,typephi);
    if ( equatn(i) == true )
        al = al + rho(i) * phival;
    else
        al = al + rho(i) * 0.5 * ( cs(i) + phival );
        % al = al + rho(i) * ( cs(i) + phival );
    end
end