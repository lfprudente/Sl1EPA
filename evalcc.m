function [c,flag] = evalcc(n,x,ind)

global problemID pairs nballs a b

% ==================================================================

if ( problemID == 1 )

    flag = 0;

    c = - x(ind);

    return

end

% ==================================================================

if ( problemID == 2 )

    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    if ( nballs == 1 )
        nck = 0;
    else
        nck = nchoosek(nballs,2);
    end
    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind; 

        c = x.r^2 - ( x.s(i) - 1 )^2 * b^2 * ( cte * x.uv{i}(1)^2 + x.uv{i}(2)^2 );

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        c = 4 * x.r^2 - a^2 * ( ( 1 + ( x.s(i) - 1 ) * cte ) * x.uv{i}(1) - ( 1 + ( x.s(j) - 1 ) * cte ) * x.uv{j}(1) )^2 ...
            - b^2 * ( x.s(i) * x.uv{i}(2) - x.s(j) * x.uv{j}(2) )^2;

    elseif ( ind >= nballs + nck + 1 && ind <= 2 * nballs + nck )

        i = ind - ( nballs + nck );

        c = - x.s(i);

    elseif ( ind >= 2 * nballs + nck + 1 && ind <= 3 * nballs + nck )

        i = ind - ( 2 * nballs + nck );

        c = x.s(i) - 1;

    else

        c = - x.r;

    end

    return
end

% ==================================================================

if ( problemID == 3 )

    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    tmp = x.A(1,1) * x.A(2,2) - x.A(2,1) * x.A(1,2);

    if ( abs(x.A(1,1)) > abs(x.A(2,2)) ) 
        center(2) = 0.5d0 * ( x.A(2,1) * x.b(1) - x.A(1,1) * x.b(2) ) / tmp;
        center(1) = - ( x.A(1,2) * center(2) + 0.5d0 * x.b(1) ) / x.A(1,1);
    else
        center(1) = 0.5d0 * ( x.A(1,2) * x.b(2) - x.A(2,2) * x.b(1) ) / tmp;
        center(2) = - ( 0.5d0 * x.b(2) + x.A(2,1) * center(1) ) / x.A(2,2);
    end
    
    if ( ind == 1 )
        c = center(1) - 10;
    elseif ( ind == 2 )
        c = - center(1) + 1;
    elseif ( ind == 3 )
        c = center(2) - 10;
    elseif ( ind == 4 )
        c = - center(2) + 1;
    end
    
    return
end