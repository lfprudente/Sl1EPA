function [dphival] = dphi(t,tau,ind)

if ( ind == 1 )
    dphival = t / sqrt( t^2 + 1/tau );
elseif ( ind == 2 )
    % a = exp(tau * t);
    % b = exp(-tau * t);
    % dphival = ( a - b ) / ( a + b );

    abst = abs(t);

    a = exp(tau*(t-abst));
    b = exp(-tau*(t+abst));

    dphival = (a-b) / (a + b);


    % coef = rho * (exp((cost_at_x-s)/epsilon)-exp((-cost_at_x-s)/epsilon))/(exp((cost_at_x-s)/epsilon)+exp((-cost_at_x-s)/epsilon));


elseif ( ind == 3 )
    if ( t >= 0.5 / tau )
        dphival = 1;
    elseif ( t <= - 0.5 / tau )
        dphival = - 1;
    else
        dphival = 2 * tau * t;
    end
elseif ( ind == 4 )
    if ( abs(t) <= 1 / tau )
        dphival = tau * t;
    else
        if ( t > 0 )
            dphival = 1;
        else
            dphival = -1;
        end
    end
elseif ( ind == 5 )
    dphival = t / sqrt( t^2 + 1/tau );
elseif ( ind == 6 )
    % a = exp(tau * t);
    % b = exp(-tau * t);
    % dphival = ( a - b ) / ( a + b );

    abst = abs(t);

    a = exp(tau*(t-abst));
    b = exp(-tau*(t+abst));

    dphival = (a-b) / (a + b);
else
    fprintf('phi: Not implemented yet.')
end