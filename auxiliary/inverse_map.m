%Source: https://github.com/mpokojovy/mvvasicek
function [R_ast, A, Sigma] = inverse_map(R, dt, solver, tol, MaxIter)
% forward_map simulates random path of a discrete multivariate Vasicek
%             process
%   [R_ast, A, Sigma] = inverse_map(R, dt, solver, tol, MaxIter)
%
% Inputs: 
%   R       : observer p-variate path of a discrete multivariate Vasicek
%             process
%   dt      : time step
%   solver  : "ml": our ML fixed point solver 
%             "ols": OLS for VAR(1)
%             "fsolve": Default Matlab solver
%             "fminsearch": Default Matlab optimizer
%             "default" or any other input: explicit solver
%   tol     : numerical solver tolerance
%   MaxIter : number of outer steps
%
% Outputs:
%   R_ast : estimated long-term mean (p x 1)
%   A     : estimated reversion speed matrix (p x p)
%   Sigma : estimated volatility matrix (p x p)
%
% Copyright (c) 2023: Michael Pokojovy, Ebenezer Nkum and
%                     Thomas M. Fullerton, Jr.

    %% Warmstart
    function [R_ast, A, Sigma] = warmstart(R, dt)
        % Initialising Parameters
        [n, p] = size(R);
        r_cur  = R(2:end, :)';       % r_cur is the r_(t_(j+1)) vector
        r_lag  = R(1:(end - 1), :)'; % r_lag is the r_(t_(j))   vector
        A1     = zeros(p, p);
        A2     = zeros(p, p);
        Sigma  = zeros(p, p);

        % Estimate for R_ast
        R_ast = mean(r_cur, 2);
       
        % Estimate for A
        for i = 1:(n - 1)
            A1 = A1 + (R_ast - r_lag(:, i))*(R_ast - r_lag(:, i))';
            A2 = A2 + (R_ast - r_lag(:, i))*(r_cur(:, i) - r_lag(:, i))';
        end
        
        A = (1/dt)*(A1\A2)';

        % Estimate for Sigma
        for i = 1:(n - 1)
            x = (r_cur(:, i) - r_lag(:, i)) - A*(R_ast - r_lag(:, i))*dt;
            Sigma = Sigma + x*x';
        end 
    
        Sigma = (1/((n - 1)*dt))*(Sigma);

        return;
    end

    %% Banach's map
    function [R_ast, A, Sigma] = banach_map(R_ast, A, Sigma, R, dt)
        % Initialising Parameters
        [n, p]  = size(R);
        r_cur   = R(2:end, :)';       % r_cur is the r_(t_(j+1)) vector
        r_lag   = R(1:(end - 1), :)'; % r_lag is the r_(t_(j))   vector
        A1      = zeros(p, p);
        A2      = zeros(p, p);
        Sigma   = zeros(p, p);

        iA = pinv(A);
        
        % Estimate for R_ast
        R_ast = mean(r_lag, 2) + iA*(R(end, :)' - R(1, :)')/(n*dt);
       
        % Estimate for A
        for i = 1:(n - 1)
            A1 = A1 + (R_ast - r_lag(:, i))*(R_ast - r_lag(:, i))';
            A2 = A2 + (r_cur(:, i) - r_lag(:, i))*(R_ast - r_lag(:, i))';
        end
        
        A = (1/dt)*A2*pinv(A1);
          
        % Estimate for Sigma
        for i = 1:(n - 1)
            x = (r_cur(:, i) - r_lag(:, i)) - A*(R_ast - r_lag(:, i))*dt;
            Sigma = Sigma + x*x';
        end 
    
        Sigma = (1/((n - 1)*dt))*(Sigma);

        return;
    end

    %% Compute log-likelihood and its gradient
    function logL = loglikelihood(theta, R, dt)
        [n, p] = size(R);

        r_cur = R(2:end, :)';   % r_cur is the r_(t_(j+1)) vector
        r_lag = R(1:end-1, :)'; % r_lag is the r_(t_(j))   vector
        
        % Parameters
        R_ast = reshape(theta(1:p), p, 1);
        A     = reshape(theta(p+1:p+p^2), p, p);
        Sigma = reshape(theta(p+p^2+1:end), p, p);
         
        % Log-likelihood
        logL = -0.5*(p*n*dt)*(log(2*pi)) - 0.5*(n*dt)*logdet(Sigma);
        
        iSigma = pinv(Sigma);

        for i = 1:(n - 1)
            x = (r_cur(:, i) - r_lag(:, i)) - A*(R_ast - r_lag(:, i))*dt;
            logL = logL - 0.5*x'*iSigma*x;
        end

        logL = -logL;

        return;
    end

    %% Computing the initial R, R_ast, A and Sigma
    [R_ast, A, Sigma] = warmstart(R, dt);

    solver = lower(solver);
    
    [n, p] = size(R);
    r_cur  = R(2:end, :)';
    r_lag  = R(1:(end - 1), :)';

    if (solver == "fsolve")
        theta0 = [R_ast; reshape(A, p*p, 1); reshape(Sigma, p*p, 1)];

        options = optimoptions('fsolve');
        options.MaxIter = MaxIter;
        options.MaxFunEvals = 1E6;
       
        sol = fsolve(@(theta) ML_est_eqs(theta, R, dt), theta0, options); % p + p^2 + p^2
         
        R_ast =  sol(1:p);
        A     =  reshape(sol((p + 1):(p + 1 + p*p - 1)), p, p);
        Sigma =  reshape(sol((p + p*p + 1):end), p, p);  
    elseif (solver == 'fminsearch')
        theta0 = [R_ast; reshape(A, p*p, 1); reshape(Sigma, p*p, 1)];

        options = optimset('MaxIter', MaxIter, 'MaxFunEvals', 1E6);
       
        sol = fminsearch(@(theta) loglikelihood(theta, R, dt), theta0, options);

        R_ast =  sol(1:p);
        A     =  reshape(sol((p + 1):(p + 1 + p*p - 1)), p, p);
        Sigma =  reshape(sol((p + p*p + 1):end), p, p); 
    elseif (solver == "ml")
        R_ast0 = R_ast;
        A0     = A;
        Sigma0 = Sigma;

        for iter = 1:MaxIter
            [R_ast, A, Sigma] = banach_map(R_ast, A, Sigma, R, dt);

            iA0     = pinv(A0);
            iSigma0 = pinv(Sigma0);

            iA     = pinv(A);
            iSigma = pinv(Sigma);

            error = mean((R_ast - R_ast0).^2, 'all') + ...
                    mean((eye(p) - iA0*A).^2, 'all') + ...
                    mean((eye(p) - A0*iA).^2, 'all') + ...
                    mean((eye(p) - iSigma0*Sigma).^2, 'all') + ...
                    mean((eye(p) - Sigma0*iSigma).^2, 'all');

            % error = mean((R_ast - R_ast0).^2, 'all') + ...
            %         mean((A - A0).^2, 'all') + ...
            %         mean((Sigma0 - Sigma).^2, 'all');   

            if (error < tol*tol)
                break
            else
                R_ast0 = R_ast;
                A0     = A;
                Sigma0 = Sigma;
            end
        end
    elseif (solver == "ols")
        model = varm(p, 1);
        est = estimate(model, R);

        A = (eye(p) - est.AR{1})/dt;
        R_ast = pinv(A)*est.Constant/dt;
        Sigma = est.Covariance/dt;
    end

    return;
end