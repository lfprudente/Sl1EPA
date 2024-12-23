function [phival] = phi(t,tau,ind)

if ( ind == 1 )
    phival = sqrt( t^2 + 1/tau );
elseif ( ind == 2 )
    % phival = log( exp(tau * t) + exp(-tau * t) ) / tau;

    abst = abs(t);

    phival = abst + log( exp( tau * ( t - abst ) ) + exp( - tau * ( t + abst ) ) ) / tau;
    
    % s = max(-cost_at_x, cost_at_x);
    % additional_cost = s + epsilon * log( exp((cost_at_x - s)/epsilon) + exp((-cost_at_x-s)/epsilon));


elseif ( ind == 3 )
    if ( t >= 0.5 / tau )
        phival = t;
    elseif ( t <= - 0.5 / tau )
        phival = - t;
    else
        phival = tau * t^2 + 0.25 / tau;
    end
elseif ( ind == 4 )
    if ( abs(t) <= 1 / tau )
        phival = 0.5 * tau * t^2;
    else
        phival = abs(t) - 0.5 / tau;
    end
elseif ( ind == 5 )
    phival = sqrt( t^2 + 1/tau ) - 1 / sqrt(tau);
elseif ( ind == 6 )
    % phival = log( exp(tau * t) + exp(-tau * t) ) / tau - log(2)/tau;

    abst = abs(t);

    phival = abst + log( exp( tau * ( t - abst ) ) + exp( - tau * ( t + abst ) ) ) / tau;
    phival = phival - log(2) / tau; 
else    
    fprintf('phi: Not implemented yet.')
end