function [rhoini] = comprhoini(m,c,f,equatn,tau,difrho)

global typephi

rhoinimin = 10^(-8);
rhoinimax = 10^8;

if ( difrho )
    for i = 1:m
        [phival] = phi(c(i),tau,typephi);
    
        if ( equatn(i) == true )
            rhoini(i) = 10 * max( 1, abs( f ) ) / max( 1, phival );
        else
            rhoini(i) = 10 * max( 1, abs( f ) ) / max( 1, 0.5 * ( c(i) + phival ) );
            % rhoini(i) = 10 * max( 1, abs( f ) ) / max( 1, c(i) + phival );
        end
    
        rhoini(i) = max( rhoinimin, min( rhoini(i), rhoinimax ) );
    end
else
    sumc = 0;
    for i = 1:m
        [phival] = phi(c(i),tau,typephi);
        
        if ( equatn(i) == true )
            sumc = sumc + phival;
        else
            sumc = sumc + 0.5 * ( c(i) + phival );
            % sumc = sumc + c(i) + phival;
        end
    end

    rhoini1 = 10 * max( 1, abs( f ) ) / max( 1, sumc );
    
    rhoini1 = max( rhoinimin, min( rhoini1, rhoinimax ) );

    rhoini(1:m) = rhoini1;

end