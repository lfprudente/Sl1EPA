clear all

global problemID pairs a b nballs typephi Points nPoints LabelPoints A

format long

problemID = 1;

% =================================================================

if ( problemID == 1 )
    
    dim_set  = [10, 20, 50, 200, 500, 1000, 2000];    % Dimension of "the Cov Matrix"
    snrset   = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0];      % Signal Strength
    deltaset = [0.1, 0.3, 0.7, 0.9];                  % Sparsity

    difrhoset = [true , false];

    for difrho = difrhoset
    
        for typephi = 1:6
                
            for dim = dim_set

                rng(2024*dim);
                
                for snr = snrset
                    
                    for delta = deltaset

                        fprintf('\n\ntypephi = %i\n',typephi)
                        fprintf('difrho  = %i\n',difrho)
                        fprintf('dim     = %i\n',dim)
                        fprintf('snr     = %f\n',snr)
                        fprintf('delta   = %f\n\n',delta)

                        % Set up data

                        T = dim;
                        samplesize = floor(delta*dim);
                        S = randsample(dim, samplesize);
                        v = zeros(dim,1);
                        v(S) = 1/sqrt(samplesize);
                        A = sqrt(snr) * v * (v.');
                        B = randn(dim)/sqrt(T);
                        for ii = 1: dim
                            B(ii,ii) = randn * 2/sqrt(T);
                        end
                        A = A+B;
                        A = (A+A')/2;      

                        % Number of variables
                        
                        n = dim;

                        % Constraints

                        m = dim;

                        E = [];
                        I = [];
                        equatn = [];
                        lambda = [];
    
                        equatn(1:m) = false;
                        lambda = zeros(m,1);
                
                        E = (equatn==true);
                        I = (equatn==false);

                        M = spherefactory(dim);
                        problem.M = M;

                        % Set initial guess

                        x0 = M.rand();
    
                        % Checking derivatives?
                        
                        checkder = false;
                        
                        % Scale the problem?
                        
                        scaling = true;
    
                        if ( checkder )
                            l(1:n) = - Inf;
                            u(1:n) =   Inf;
                            checkd(n,m,x0,l,u)
                        end
            
                        % Set the feseabilty, optimality, and complementarity tolerances
                        
                        epsfeas   = 10^(-4);
                        epsopt    = 10^(-4);
                        epscompl  = 10^(-4);

                        % Call the solver
                                
                       [x,lambda,f,outiter,initer,csupn,snorm,nlmannorm,maxrho,tau,time,nevalal,nevalnal,alinfo] = Sl1EPA(n,m,x0,equatn,lambda,difrho,scaling,problem,epsopt,epsfeas,epscompl);

                    end
                end
            end
        end
    end
end

% ==================================================================

if ( problemID == 2 )

    balls = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];

    difrho = false;

    for j = 1:length(balls)

        % Set the number of balls to be packed

        nballs = balls(j);         

        for typephi = 1:6

            rng(2024*nballs);
   
            % Print solution?
            
            print = true;
        
            % Problem data
                
            a = 2;
            b = 1;
        
            % Number of variables
        
            n = 3 * nballs + 1;
        
            % Constraints
        
            if ( nballs == 1 )
                m = 3 * nballs + 1;
            else
                m = nchoosek(nballs,2) + 3 * nballs + 1;
            end
    
            equatn = [];
            lambda = [];
            E = [];
            I = [];
    
            equatn(1:m) = false;
            lambda = zeros(m,1);
            
            E = (equatn==true);
            I = (equatn==false);
        
            sphere = spherefactory(2);
            manifold = productmanifold(struct('uv',powermanifold(sphere, nballs),'s', euclideanfactory(nballs),'r', euclideanfactory(1)));
            problem.M = manifold;
    
            if ( nballs == 1 )
                pairs = 0;
            else
                pairs = nchoosek(1:nballs, 2);
            end
        
            % Checking derivatives?
            
            checkder = false;
            
            % Scale the problem?
            
            scaling = true;
        
            % Set initial guess
    
            x = [];
        
            for i = 1:nballs
                x0 = 2 * rand(2,1) - 1;
                x0 = x0 / norm(x0);
                x.uv{i}(:,1) = x0;
            end
            x.uv = x.uv';
            x.s = rand(nballs,1);
            x.r = rand;
            
            if ( checkder )
                xa = cell2mat(x.uv);
                xb = [];
                xb(1:nballs,1) = xa(1:2:end,1);
                xb(nballs+1:2*nballs,1) = xa(2:2:end,1);
                xz = [xb; x.s; x.r];
                checkd(n,m,xz,-Inf(n,1),Inf(n,1))
            end               
            x0 = x;    

            % Set the feseabilty, optimality, and complementarity tolerances
                        
            epsfeas   = 10^(-4);
            epsopt    = 10^(-4);
            epscompl  = 10^(-4);

            % Call the solver

            fprintf('\n\nnballs = %i \n',nballs)
            fprintf('Typephi = %i \n\n', typephi)
        
            [x,lambda,f,outiter,initer,csupn,snorm,nlmannorm,maxrho,tau,time,nevalal,nevalnal,alinfo] = Sl1EPA(n,m,x0,equatn,lambda,difrho,scaling,problem,epsopt,epsfeas,epscompl);
            
             if ( print ) printfigure(x), end;
   
        end
    end
end


% ==================================================================

if ( problemID == 3 )

    difrho = false;

    for Ptype = 1:4

        for typephi = 1:6

            rng(2024*Ptype);
        
            % Number of variables
            
            n = 6;
            
            % Constraints  

            m = 4;
            
            equatn(1:m) = false;
            lambda = zeros(m,1);
            
            E = (equatn==true);
            I = (equatn==false);
        
            manifold = productmanifold(struct('A',sympositivedefinitefactory(2),'b', euclideanfactory(2)));;
            problem.M = manifold;
            
            % Checking derivatives?
            
            checkder = false;
            
            % Scale the problem?
            
            scaling = true;
        
            % Print solution?

            print = false;
        
            % Set the problem data

            nPoints = 10000;
            Points  = -10 + 20 * rand(2,nPoints);
        
            for i = 1:nPoints
                if (  ( Ptype == 1 && norm(Points(:,i)) <= 7 ) || ...
                      ( Ptype == 2 && norm(Points(:,i),Inf) <= 7 ) || ...
                      ( Ptype == 3 && abs(Points(1,i)) <= 7 && abs(Points(2,i)) <= 3.5 ) || ...
                      ( Ptype == 4 && 2 * Points(1,i) - Points(2,i) <= 7 && - Points(1,i) + 2 * Points(2,i) <= 7 && - Points(1,i) - Points(2,i) <= 7 ) )
                    LabelPoints(i) = 1;
                else
                    LabelPoints(i) = 0;
                end
            end
        
            % Set initial guess
        
            x.A = rand * eye(2);
            x.b = rand(2,1);
      
            if ( checkder )
                xz = [x.A(:);x.b(:)];
                l(1:n) = - Inf;
                u(1:n) =   Inf;
                checkd(n,m,xz,l,u)
            end
        
            % Set the feseabilty, optimality, and complementarity tolerances
                        
            epsfeas   = 10^(-6);
            epsopt    = 10^(-6);
            epscompl  = 10^(-6);
    
            % Call the solver

            fprintf('\n\nPtype   = %i \n',Ptype)
            fprintf('Typephi = %i \n\n', typephi)
  
            [x,lambda,f,outiter,initer,csupn,snorm,nlmannorm,maxrho,tau,time,nevalal,nevalnal,alinfo] = Sl1EPA(n,m,x,equatn,lambda,difrho,scaling,problem,epsopt,epsfeas,epscompl);
            
            if ( print ) printfigure(x); end
        end

    end
end