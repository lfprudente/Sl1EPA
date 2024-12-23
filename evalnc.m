function [nc,flag] = evalnc(n,x,ind)

global problemID pairs nballs a b

% ==================================================================

if ( problemID == 1 )

    flag = 0;

    nc = zeros(n,1);

    nc(ind) = - 1;

    return

end

% ==================================================================

if ( problemID == 2 )

    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    nc = zeros(n, 1);

    if ( nballs == 1 )
        nck = 0;
    else
        nck = nchoosek(nballs,2);
    end


    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind;

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Pr  = 3 * nballs + 1;

        nc(Pui) = - ( x.s(i) - 1 )^2 * b^2 * cte * 2 * x.uv{i}(1);
        nc(Pvi) = - ( x.s(i) - 1 )^2 * b^2 * 2 * x.uv{i}(2);
        nc(Psi) = - 2 * ( x.s(i) - 1 )* b^2 * ( cte * x.uv{i}(1)^2 + x.uv{i}(2)^2 );
        nc(Pr)  = 2 * x.r;

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Puj = j;
        Pvj = j + nballs;
        Psj = j + 2 * nballs;
        Pr  = 3 * nballs + 1;

        xxi = ( 1 + ( x.s(i) - 1 ) * cte ) * x.uv{i}(1);
        xxj = ( 1 + ( x.s(j) - 1 ) * cte ) * x.uv{j}(1);
        yyi = x.s(i) * x.uv{i}(2);
        yyj = x.s(j) * x.uv{j}(2);

        nc(Pui) = - 2 * a^2 * ( xxi - xxj ) * ( 1 + ( x.s(i) - 1 ) * cte );
        nc(Pvi) = - 2 * b^2 * ( yyi - yyj ) * x.s(i);
        nc(Psi) = - 2 * a^2 * ( xxi - xxj ) * cte * x.uv{i}(1) - 2 * b^2 * ( yyi - yyj ) * x.uv{i}(2);
        nc(Puj) =   2 * a^2 * ( xxi - xxj ) * ( 1 + ( x.s(j) - 1 ) * cte );
        nc(Pvj) =   2 * b^2 * ( yyi - yyj ) * x.s(j);
        nc(Psj) =   2 * a^2 * ( xxi - xxj ) * cte * x.uv{j}(1) + 2 * b^2 * ( yyi - yyj ) * x.uv{j}(2);
        nc(Pr)  =   8 * x.r;

    elseif ( ind >= nballs + nck + 1 && ind <= 2 * nballs + nck )

        i   = ind - ( nballs + nck );
        Psi = i + 2 * nballs;

        nc(Psi) = - 1;

    elseif ( ind >= 2 * nballs + nck + 1 && ind <= 3 * nballs + nck )

        i   = ind - ( 2 * nballs + nck );
        Psi = i + 2 * nballs;

        nc(Psi) = 1;
    else
        Pr  = 3 * nballs + 1;

        nc(Pr) = - 1;
    end

    return
end

% ==================================================================


if ( problemID == 3 )

    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    % Initialize the gradient vector with zeros
    nc = zeros(n, 1);

    tmp = x.A(1,1) * x.A(2,2) - x.A(2,1) * x.A(1,2);
    dtmpdA11 = x.A(2,2);
    dtmpdA22 = x.A(1,1);
    dtmpdA21 = - x.A(1,2);
    dtmpdA12 = - x.A(2,1);
    dtmpdb1  = 0;
    dtmpdb2  = 0;

    if ( abs(x.A(1,1)) > abs(x.A(2,2)) ) 
        num2 = 0.5d0 * ( x.A(2,1) * x.b(1) - x.A(1,1) * x.b(2) );
        dnum2dA11 = - 0.5d0 * x.b(2);
        dnum2dA12 = 0;
        dnum2dA21 = 0.5d0 * x.b(1);
        dnum2dA22 = 0;
        dnum2db1  = 0.5d0 * x.A(2,1);
        dnum2db2  = - 0.5d0 * x.A(1,1);

        center(2) = num2 / tmp;

        dcenter2dA11 = dnum2dA11 / tmp - num2 * dtmpdA11 / tmp^2;
        dcenter2dA12 = dnum2dA12 / tmp - num2 * dtmpdA12 / tmp^2;
        dcenter2dA21 = dnum2dA21 / tmp - num2 * dtmpdA21 / tmp^2;
        dcenter2dA22 = dnum2dA22 / tmp - num2 * dtmpdA22 / tmp^2;
        dcenter2db1 = dnum2db1 / tmp - num2 * dtmpdb1 / tmp^2;
        dcenter2db2 = dnum2db2 / tmp - num2 * dtmpdb2 / tmp^2;

        num1 = - ( x.A(1,2) * center(2) + 0.5d0 * x.b(1) );

        dnum1dA11 = - x.A(1,2) * dcenter2dA11;
        dnum1dA12 = - center(2) - x.A(1,2) * dcenter2dA12;
        dnum1dA21 = - x.A(1,2) * dcenter2dA21;
        dnum1dA22 = - x.A(1,2) * dcenter2dA22;
        dnum1db1  = - x.A(1,2) * dcenter2db1 - 0.5d0;
        dnum1db2  = - x.A(1,2) * dcenter2db2;

        center(1) = num1 / x.A(1,1);

        dcenter1dA11 = dnum1dA11 / x.A(1,1) - num1 / x.A(1,1)^2;
        dcenter1dA12 = dnum1dA12 / x.A(1,1);
        dcenter1dA21 = dnum1dA21 / x.A(1,1);
        dcenter1dA22 = dnum1dA22 / x.A(1,1);
        dcenter1db1  = dnum1db1 / x.A(1,1);
        dcenter1db2  = dnum1db2 / x.A(1,1);
    else
        num1 = 0.5d0 * ( x.A(1,2) * x.b(2) - x.A(2,2) * x.b(1) );

        dnum1dA11 = 0;
        dnum1dA12 = 0.5d0 * x.b(2);
        dnum1dA21 = 0;
        dnum1dA22 = - 0.5d0 * x.b(1);
        dnum1db1  = - 0.5d0 * x.A(2,2);
        dnum1db2  = 0.5d0 * x.A(1,2);

        center(1) = num1 / tmp;

        dcenter1dA11 = dnum1dA11 / tmp - num1 * dtmpdA11 / tmp^2;
        dcenter1dA12 = dnum1dA12 / tmp - num1 * dtmpdA12 / tmp^2;
        dcenter1dA21 = dnum1dA21 / tmp - num1 * dtmpdA21 / tmp^2;
        dcenter1dA22 = dnum1dA22 / tmp - num1 * dtmpdA22 / tmp^2;
        dcenter1db1  = dnum1db1 / tmp - num1 * dtmpdb1 / tmp^2;
        dcenter1db2  = dnum1db2 / tmp - num1 * dtmpdb2 / tmp^2;

        num2 = - ( 0.5d0 * x.b(2) + x.A(2,1) * center(1) );

        dnum2dA11 = - x.A(2,1) * dcenter1dA11;
        dnum2dA12 = - x.A(2,1) * dcenter1dA12;
        dnum2dA21 = - center(1) - x.A(2,1) * dcenter1dA21;
        dnum2dA22 = - x.A(2,1) * dcenter1dA22;
        dnum2db1  = - x.A(2,1) * dcenter1db1;
        dnum2db2  = - x.A(2,1) * dcenter1db2 - 0.5d0;


        center(2) = num2 / x.A(2,2);

        dcenter2dA11 = dnum2dA11 / x.A(2,2);
        dcenter2dA12 = dnum2dA12 / x.A(2,2);
        dcenter2dA21 = dnum2dA21 / x.A(2,2);
        dcenter2dA22 = dnum2dA22 / x.A(2,2) - num2 / x.A(2,2)^2;
        dcenter2db1  = dnum2db1 / x.A(2,2);
        dcenter2db2  = dnum2db2 / x.A(2,2);

    end
    
    if ( ind == 1 )
        nc(1) = dcenter1dA11;
        nc(2) = dcenter1dA21;
        nc(3) = dcenter1dA12;
        nc(4) = dcenter1dA22;
        nc(5) = dcenter1db1;
        nc(6) = dcenter1db2;
    elseif ( ind == 2 )
        nc(1) = dcenter1dA11;
        nc(2) = dcenter1dA21;
        nc(3) = dcenter1dA12;
        nc(4) = dcenter1dA22;
        nc(5) = dcenter1db1;
        nc(6) = dcenter1db2;
        
        nc = - nc;
    elseif ( ind == 3 )
        nc(1) = dcenter2dA11;
        nc(2) = dcenter2dA21;
        nc(3) = dcenter2dA12;
        nc(4) = dcenter2dA22;
        nc(5) = dcenter2db1;
        nc(6) = dcenter2db2;
    elseif ( ind == 4 )
        nc(1) = dcenter2dA11;
        nc(2) = dcenter2dA21;
        nc(3) = dcenter2dA12;
        nc(4) = dcenter2dA22;
        nc(5) = dcenter2db1;
        nc(6) = dcenter2db2;

        nc = - nc;
    end
end