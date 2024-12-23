function [x,lambda,f,outiter,initer,csupn,snorm,nlmannorm,maxrho,tau,time,nevalal,nevalnal,alinfo] = Sl1EPA(n,m,x,equatn,lambda,difrho,scaling,problem,epsopt,epsfeas,epscompl)

global evalfun evalgrad nevalal nevalnal typephi

%  Parameters

epsfeas12 = sqrt( epsfeas );
epsopt12  = sqrt( epsopt );

rhomax  = 10^20;
rhofrac = 0.5;
rhomult = 10;
taufrac = 10;
tau = 1;

maxoutit = 50;
maxinnit = 1000;

E = ( equatn == true  );
I = ( equatn == false );

if ( scaling == true )
    [g,~] = evalg(n,x);
    sc.f = 1 / max( 1, norm(g,Inf) );
    
    for i = 1:m
        [nc,~] = evalnc(n,x,i);
        sc.c(i) = 1 / max( 1, norm(nc,Inf) );
    end
else
    sc.f = 1;
    sc.c(1:m) = 1;
end

% ==================================================================
% Print initial information
% ==================================================================

fprintf('----------------------------------------------------------------------\n')
fprintf('Smooth l1-exact penalty algorithm for Riemannian optimization problems\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('Number of variables                : %i \n',n)
fprintf('Number of constraints              : %i \n\n',m)
fprintf('Optimality tolerance               : %.0e \n',epsopt)
fprintf('Feasibility tolerance              : %.0e \n',epsfeas)
fprintf('Complementarity tolerance          : %.0e \n\n',epscompl)
fprintf('Smoothing function                 : %i \n',typephi)
if ( scaling == true )
    fprintf('Objective function scale factor    : %.0e \n',sc.f)
    fprintf('Smallest constraints scale factor  : %.0e \n',min(sc.c))
end

% Start timing

tic;

%  ==================================================================
%  Initialization
%  ==================================================================

%  Counters

outiter  = 0;
initer   = 0;
evalfun  = 0;
evalgrad = 0;

nevalal  = 0;
nevalnal = 0;

ISerror  = 0;

% Compute objective function value and gradient

[f,fs,~] = sevalf(n,x,sc,scaling);
evalfun = evalfun + 1;

[~,gs,~] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

% Compute constraints 

c  = [];
cs = [];
for i = 1:m
    [c(i),cs(i),~] = sevalc(n,x,i,sc,scaling);
    evalfun = evalfun + 1;
end
c  = c';
cs = cs';

% Compute complementarity and feasibility violations

csupn  = max( norm(abs(c(E)),Inf), norm(max(c(I),0),Inf) );
cssupn = max( norm(abs(cs(E)),Inf), norm(max(cs(I),0),Inf) );
snorm  = norm(min(-cs(I),lambda(I)),Inf);

% Compute the Euclidian Lagrangian gradient

[nl,ni] = evalnlni(n,m,x,lambda,lambda,gs,sc,scaling);
evalgrad = evalgrad + m;

% Convert the Euclidean Lagrangian gradient to the Riemannian Lagrangian gradient

nlman = problem.M.egrad2rgrad(x,nl);

% Compute the norm of the Riemannian Lagrangian gradient

nlmannorm = problem.M.norm(x,nlman);

% Convert the Euclidean squared-infeasibility gradient to the Riemannian squared-infeasibility gradient

niman = problem.M.egrad2rgrad(x,ni);

% Compute squared-infeasibility gradient norm

nimannorm = problem.M.norm(x,niman);

% Save best solution

fb       = f;
csupnb   = csupn;
snormb   = snorm;
nlmannormb = nlmannorm;
xb       = x;
lambdab  = lambda;

% Option structures from Manopt

warning('off', 'manopt:getHessian:approx');
options.verbosity = 0;

% ==================================================================
% Main loop
% ==================================================================

while (1)

    % ==================================================================
    % Print information of this iteration
    % ==================================================================
    
    % Print information    

    if ( scaling == true )
        if ( mod(outiter,10) == 0 )
            fprintf('\n')
            fprintf('%3s %-6s   %-6s  %-9s %-6s  %-6s  %+9s %-6s %-4s   %-6s %-5s \n','out','penalt','smooth','objective','infeas','scaled','scaled','comple','norm','|Grad|','inner')
            fprintf('%3s %-6s   %-6s  %-9s %-6s  %-6s  %-6s %-6s %-4s %-6s %-5s \n','ite','param ','param ','function ','ibilty','obj-funct','infeas','mentar','graLag','infeas','stop ')
        end
        if ( outiter == 0 )
            fprintf('%3i   -     %6.1e %10.3e %6.0e %10.3e  %6.0e %6.0e %6.0e %6.0e  - \n',outiter,tau,f,csupn,fs,cssupn,snorm,nlmannorm,nimannorm)
        else
            fprintf('%3i %6.0e  %6.1e %10.3e %6.0e %10.3e  %6.0e %6.0e %6.0e %6.0e  %i \n',outiter,maxrho,tau,f,csupn,fs,cssupn,snorm,nlmannorm,nimannorm,flagIS)
        end

    else
        if ( mod(outiter,10) == 0 )
            fprintf('\n')
            fprintf('%3s %-6s  %-6s %-9s %-6s %-6s %-4s   %-6s %-5s \n','out','penalt','smooth','objective','infeas','comple','norm','|Grad|','inner')
            fprintf('%3s %-6s  %-6s %-9s %-6s %-6s %-4s %-6s %-5s \n','ite','param ','param ','function ','ibilty','mentar','graLag','infeas','stop ')
        end
        if ( outiter == 0 )
            fprintf('%3i   -    %6.1e %10.3e %6.0e %6.0e %6.0e %6.0e - \n',outiter,tau,f,csupn,snorm,nlmannorm,nimannorm)
        else
            fprintf('%3i %6.0e  %6.1e %10.3e %6.0e %6.0e %6.0e %6.0e %i \n',outiter,maxrho,tau,f,csupn,snorm,nlmannorm,nimannorm,flagIS)
        end
    end
    
    % ==================================================================
    % Test stopping criteria
    % ==================================================================
    
    % Test feasibility, optimality and complementarity
    
    if ( csupn <= epsfeas &&  snorm <= epscompl && nlmannorm <= epsopt )
     
        alinfo = 0;
        
        % Stop timing 
     
        time = toc;
    
        % Print information
    
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)        
 
        return
    end
    
    % Test whether we are at an infeasible point that is stationary for
    % the sum of the squared infeasibilities
    
%     if ( csupn > efstain && nimannorm <= eostain ) 
%     
%         alinfo = 1;
%         
%         % Stop timing 
%         
%         time = toc;
%         
%         % Print information
%         
%         fprintf('\n')
%         fprintf('It seems that a stationary-of-the-infeasibility probably infeasible point was found.\n')
%         fprintf('Whether the final iterate is a solution or not requires further analysis.\n')
%         fprintf('Number of function evaluations: %i\n',evalfun)
%         fprintf('Number of gradient evaluations: %i\n',evalgrad)
%         fprintf('CPU time(s)                   : %.1f \n',time)
%         
%         return
%     end
    
    % Test whether the penalty parameter is too large
    
    if ( outiter > 0 && maxrho > rhomax )
    
        alinfo = 2;
        
        % Stop timing 
        
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('The penalty parameter is too large. The problem may be \n')
        fprintf('infeasible or badly scaled. Further analysis is required.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)
        
        return
    end
    
    % Test whether the number of iterations is exhausted
    
    if ( outiter >= maxoutit )
        alinfo = 3;
        
        % Stop timing 
        
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Maximum of iterations reached. The feasibility-complementarity and \n')
        fprintf('optimality tolerances could not be achieved.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)
        
        return
    end

    % Test for consecutive fails in the Inner Solver

    if ( outiter > 0 && ISerror >= 2 &&  csupn <= epsfeas )
        alinfo = 4;
        
        % Stop timing 
        
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('The subproblem can not be solved. The complementarity and \n')
        fprintf('optimality tolerances could not be achieved.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)
        
        return
    end
    
    % ==================================================================
    % Iteration
    % ==================================================================
    
    outiter = outiter + 1;
    
    % ==================================================================
    % Set penalty parameter
    % ==================================================================
    
    if ( outiter == 1 ) 

        [rho] = comprhoini(m,cs,fs,equatn,tau,difrho);

    else
        if ( difrho )
            for i = 1:m
                if  ( ( equatn(i) == true  && abs(cs(i))   > rhofrac * abs(csprev(i))   ) || ...
                      ( equatn(i) == false && max(0,cs(i)) > rhofrac * max(0,csprev(i)) ) ) 
                    rho(i) = rhomult * rho(i);
                end
            end
        else
            if ( max( norm(cs(E),inf),norm(max(0,cs(I)),inf) ) > rhofrac * max( norm(csprev(E),inf),norm(max(0,csprev(I)),inf) ) )
                rho = rhomult * rho;
            end
        end

    end

    maxrho = max(rho);
    
    % ==================================================================
    % Solve the augmented Lagrangian subproblem
    % ==================================================================
    
    % Set optimality requeriment for the subproblem
    
    if ( outiter == 1 ) 
        epsopk = sqrt( epsopt );
    elseif ( csupn <= epsfeas12 && nlmannorm <= epsopt12 )
        epsopk = min( rhofrac * nlmannorm, 0.1 * epsopk );
        epsopk = max( epsopk, epsopt );
    end

    % Problem and option structures from Manopt

    problem.cost  = @(x,store) evalal(x,store,n,m,equatn,tau,rho,sc,scaling);
    problem.egrad = @(x,store) evalnal(x,store,n,m,equatn,tau,rho,sc,scaling);

    options.tolgradnorm = epsopk;
    
    if ( outiter == 1 )
        options.maxiter = 100;
    else
        options.maxiter = maxinnit;
    end
    options.maxiter = maxinnit;

    if ( m == 0 )
        options.tolgradnorm = epsopt;
        options.maxiter = 100000;
    end

    % checkgradient(problem); pause;

    [x,~,info,options] = rlbfgs(problem,x,options);

    if ( info(end).gradnorm <= epsopk || outiter == 1 )
        flagIS = 0;
    elseif ( info(end).iter >= options.maxiter )
        flagIS = 1;
    else
        flagIS = 2;
    end

    if ( flagIS ~= 0 )
        ISerror = ISerror + 1;
    else
        ISerror = 0;
    end

    initer = initer + info(end).iter;

    % ==================================================================
    % Prepare for the next iteration
    % ==================================================================
    
    % Compute objective function value and gradient
    
    [f,fs,~] = sevalf(n,x,sc,scaling);
    evalfun = evalfun + 1;

    [~,gs,~] = sevalg(n,x,sc,scaling);
    evalgrad = evalgrad + 1;
    
    % Compute constraints 

    csprev = cs;
    
    c  = [];
    cs = [];
    for i = 1:m
        [c(i),cs(i),~] = sevalc(n,x,i,sc,scaling);
        evalfun = evalfun + 1;
    end
    c  = c';
    cs = cs';
    
    % Compute feasibility violation
    
    csupn  = max( norm(abs(c(E)),Inf), norm(max(c(I),0),Inf) );
    cssupn = max( norm(abs(cs(E)),Inf), norm(max(cs(I),0),Inf) );
    
    % Update Lagrange multipliers approximation

    for i = 1:m
        
        [dphival] = dphi(cs(i),tau,typephi);

        rhoesc = rho / max(rho);

        if ( equatn(i) == true )
            lambda(i)    = rho(i) * dphival;
            lambdaesc(i) = rhoesc(i) * dphival;
        else
            lambda(i)    = rho(i) * 0.5 * ( 1 + dphival );
            lambdaesc(i) = rhoesc(i) * 0.5 * ( 1 + dphival );
        end

    end
    
    % Compute complementarity violation    
    
    snorm = norm(min(-cs(I),lambda(I)),Inf);

    % Compute the Euclidian Lagrangian gradient
    
    [nl,ni] = evalnlni(n,m,x,lambda,lambdaesc,gs,sc,scaling);
    evalgrad = evalgrad + m;
    
    % Convert the Euclidean Lagrangian gradient to the Riemannian Lagrangian gradient
    
    nlman = problem.M.egrad2rgrad(x,nl);
    
    % Compute the norm of the Riemannian Lagrangian gradient
    
    nlmannorm = problem.M.norm(x,nlman);

    % Convert the Euclidean squared-infeasibility gradient to the Riemannian squared-infeasibility gradient

    niman = problem.M.egrad2rgrad(x,ni);
    
    % Compute squared-infeasibility gradient norm
    
    nimannorm = problem.M.norm(x,niman);

    % Update smoothing parameter

    tau = taufrac * tau;
    
    % Save best solution
    
    if ( ( csupnb > epsfeas && csupn < csupnb ) || ...
     ( csupnb <= epsfeas && csupn <= epsfeas && f < fb ) ) 
    
        fb       = f;
        csupnb   = csupn;
        snormb   = snorm;
        nlmannormb = nlmannorm;
        xb       = x;
        lambdab  = lambda;
    
    end

    % ==================================================================
    % Iterate
    % ==================================================================

end

% ==================================================================
% End of main loop
% ==================================================================