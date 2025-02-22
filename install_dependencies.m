function install_dependencies()

    % Installing dependencies

    % Define path
    manoptDir = fullfile(pwd, 'manopt');
    
    % Run importmanopt.m to set up Manopt
    fprintf('\n\n=================================\n');
    fprintf('Installing Manopt...\n');

    try
        currentDir = pwd;
        cd(manoptDir);
        importmanopt();
        cd(currentDir);
    catch
        fprintf('\nFailed to install Manopt.\n');
        return
    end

    % Test Manopt installation
    test_manopt = input('\nDo you want to test Manopt installation with a sample problem? [Y/N] ', 's');
    if strcmpi(test_manopt, 'Y') || strcmpi(test_manopt, 'y')
        fprintf('Running sample problem for Manopt...\n');
        [flagManopt] = run_sample_problem_manopt();

        if ( flagManopt == 0 )
            fprintf('\nManopt converged successfully\n');
        else
            fprintf('\nManopt did not converge as it should.\n');
        end
    end
    
    % Save the MATLAB path
    response = input('Save paths for future Matlab sessions? [Y/N] ', 's');
    if strcmpi(response, 'Y')
        failed = savepath();
        if ~failed
            fprintf('Path saved successfully: no need to call install_dependencies next time.\n');
        else
            fprintf(['Failed to save the path. You might need administrator privileges \nto write on pathdef.m?\nPath not saved: ' ...
                     'please re-call install_dependencies next time.\n']);
        end
    else
        fprintf('Path not saved: please re-call install_dependencies next time.\n');
    end
    
    % Confirm installation
    fprintf('\n\nDependencies installed successfully.\n');
end

function [status] = run_sample_problem_ASA()
    n = 100;
    lo = zeros(n,1);    % lower bound
    hi = ones(n,1);     % upper bound
    x  = ones(n,1);     % initial guess; necessary so that it knows the size of the problem
    tol = 10^(-8);      % tolerance to optimality

    fcn  = 'asa_fcn';
    grad = 'asa_grad'; 
    fcnGrad = 'asa_fcnGrad'; % optional

    opts = [];
    opts.PrintParms = false;
    opts.PrintLevel = 1;
    opts.PrintFinal = true;
    CGopts = struct('PrintParms',false);

    fgopt.n = n;
    
    [x,status,statistics] = asa_wrapper( x, lo, hi, tol,fcn, grad, fcnGrad, opts, CGopts, fgopt);
    
    fprintf('\n');
    fprintf('Sample problem for ASA completed.\n');
end

function [status] = run_sample_problem_manopt()
    % Generate random problem data.

    rng(2024);
    n = 1000;
    A = randn(n);
    A = .5*(A+A');

    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;

    x = manifold.rand();

    % Define the problem cost function and its Euclidean gradient.
    problem.cost  = @(x) -x'*(A*x);
    problem.egrad = @(x) -2*A*x;      % notice the 'e' in 'egrad' for Euclidean

    options.tolgradnorm = 10^(-4);

    % Numerically check gradient consistency (just once, optional).
    % checkgradient(problem); pause;

    % Solve.
    [x, xcost, info, options] = rlbfgs(problem,x,options);
    
    status = - 1;
    if ( info(end).gradnorm <= options.tolgradnorm )
        status = 0;
    end

    % Display some statistics.
%     figure;
%     semilogy([info.iter], [info.gradnorm], '.-');
%     xlabel('Iteration number');
%     ylabel('Norm of the Riemannian gradient of f');
    
    fprintf('\n');
    fprintf('Sample problem for Manopt completed.\n');
end
