function [R_ast, A, Sigma] = inverse_map_rob(R, dt, OutlierFraction, solver, tol, MaxIter)
% forward_map simulates random path of a discrete multivariate Vasicek
%             process
%   [R_ast, A, Sigma] = inverse_map(R, dt, solver, tol, MaxIter)
%
% Inputs: 
%   R       : observed p-variate path of a discrete multivariate Vasicek
%             process
%   dt      : time step
%   solver  : "mtle": Minimum Trimmed Likelihood Estimation
%             "mlts": MLTS for VAR(1)
%   tol     : numerical solver tolerance
%   MaxIter : number of outer steps
%
% Outputs:
%   R_ast : estimated long-term mean (p x 1)
%   A     : estimated reversion speed matrix (p x p)
%   Sigma : estimated volatility matrix (p x p)
%
% Copyright (C) 2024: Andrews T. Anum, Michael Pokojovy, Ebenezer Nkum and
%                     Thomas M. Fullerton, Jr.

    [n, p] = size(R);

    if (nargin < 6)
        MaxIter = 1000;
    elseif (nargin < 5)
        tol = 1E-6;
    elseif (nargin < 4)
        solver = 'mlts';
    elseif (nargin < 3)
        OutlierFraction = ((n - 1) - floor(((n - 1) + p + 1)/2))/n;
    end

    %% Banach's map
    function [R_ast, A, Sigma] = banach_map(R_ast, A, Sigma)
        %iA     = pinv(A);

        % [U,S,V] = svd(A);
        % %d_mtle = pmin(pmax(A, max([min(d_mlts), 1E-3*max(d_mlts), 1E-3*max(diag(S))])), max(d_mlts));
        % d_mtle = max(diag(S), max([0.5*min(d_mlts), 1E-3*max(d_mlts), 1E-3*max(diag(S))]));
        % iA = U*diag(1./d_mtle)*V';

        [U,S,V] = svd(A);
        d_mtle = diag(S);
        d_mtle = min(d_mtle, 2*max(d_mlts));
        d_mtle = max(d_mtle, max([0.5*min(d_mlts), 1E-3*max(d_mlts), 1E-3*max(d_mtle)]));
        % iA = U*diag(1./d_mtle)*V';
        iA = V*diag(1./d_mtle)*U';

        iSigma = pinv(Sigma);

        A1    = zeros(p, p);
        A2    = zeros(p, p);
        Sigma = zeros(p, p);
        mah   = zeros(n - 1, 1);

        for i = 1:(n - 1)
            x = (r_cur(:, i) - r_lag(:, i)) - A*(R_ast - r_lag(:, i))*dt;
            mah(i) = x'*iSigma*x;
        end

        [~, I] = sort(mah, 'ascend');
        I = reshape(I, 1, n - 1);
        I = I(1:h);

        % Estimate for R_ast
        %R_ast = mean(r_lag(:, I), 2) + iA*(R(end, :)' - R(1, :)')/(h*dt);
        R_ast = zeros(p, 1);
        for i = I
            R_ast = R_ast + (r_cur(:, i) - r_lag(:, i));
        end
        R_ast = mean(r_lag(:, I), 2) + iA*R_ast*(1/(h*dt));  %/(h*dt);

        % Estimate for A
        for i = I
            A1 = A1 + (R_ast - r_lag(:, i))*(R_ast - r_lag(:, i))';
            A2 = A2 + (r_cur(:, i) - r_lag(:, i))*(R_ast - r_lag(:, i))';
        end
        
        A = (1/dt)*A2*pinv(A1);

        [U,S,V] = svd(A);
        d_mtle = diag(S);
        d_mtle = min(d_mtle, 2*max(d_mlts));
        d_mtle = max(d_mtle, max([0.5*min(d_mlts), 1E-3*max(d_mlts), 1E-3*max(d_mtle)]));
        A = U*diag(d_mtle)*V';
          
        % Estimate for Sigma
        for i = I
            x = (r_cur(:, i) - r_lag(:, i)) - A*(R_ast - r_lag(:, i))*dt;
            Sigma = Sigma + x*x';
        end 
    
        Sigma = (1/((h - p + 1)*dt))*(Sigma);

        return;
    end

    solver = lower(solver);
    
    h = floor(n*(1 - OutlierFraction)) + 1;

    r_cur  = R(2:end, :)';
    r_lag  = R(1:(end - 1), :)';

    %% Compute MLTS estimates
    x = R(1:end - 1, :);
    y = R(2:end, :);
    x = [ones(n - 1, 1) x];

    [beta, sigma, dres] = mlts(x, y, OutlierFraction);

    [~, I] = sort(dres, 'ascend');
    I = reshape(I, 1, n - 1);
    I = I(1:h);

    A = (eye(p) - beta(2:end, :)')/dt;
    R_ast = pinv(A)*beta(1, :)'/dt;
    %R_ast = mean(r_lag(:, I), 2);
    Sigma = sigma/dt;

    [~,S,~] = svd(A);
    d_mlts = diag(S);
    

    % A1    = zeros(p, p);
    % A2    = zeros(p, p);
    % 
    % for i = I
    %     A1 = A1 + (R_ast - r_lag(:, i))*(R_ast - r_lag(:, i))';
    %     A2 = A2 + (r_cur(:, i) - r_lag(:, i))*(R_ast - r_lag(:, i))';
    % end
    % A = (1/dt)*A2*pinv(A1);

    %% If selected, apply MTLE using MLTS as warmstart
    if (solver == "mlts")
        [beta, sigma] = mlts(x, y, OutlierFraction);

        A = (eye(p) - beta(2:end, :)')/dt;
        R_ast = pinv(A)*beta(1, :)'/dt;
        Sigma = sigma/dt;
    else
        R_ast0 = R_ast;
        A0     = A;
        Sigma0 = Sigma;

        for iter = 1:MaxIter
            [R_ast, A, Sigma] = banach_map(R_ast, A, Sigma);

            iSigma = pinv(0.5*(Sigma0 + Sigma));

            err_R = sqrt((1/p)*(R_ast - R_ast0)'*iSigma*(R_ast - R_ast0));
            err_A = sqrt(mean((A - A0).^2, 'all'));

            iSigma0 = pinv(Sigma0);
            iSigma  = pinv(Sigma);

            err_S = sqrt(mean((eye(p) - iSigma0*Sigma).^2, 'all') + ...
                         mean((eye(p) - Sigma0*iSigma).^2, 'all'));

            error = err_R/(1.0 + err_R) + err_A/(1.0 + err_A) + ...
                    err_S/(1.0 + err_S);

            if (error < tol)
                break
            else
                R_ast0 = R_ast;
                A0     = A;
                Sigma0 = Sigma;
            end
        end

        Sigma = (1 - OutlierFraction)/...
                chi2cdf(chi2inv(1 - OutlierFraction, p), p + 2)*Sigma;

        iA    = pinv(A, 0.1*norm(A));
        R_ast = zeros(p, 1);

        for i = I
            R_ast = R_ast + (r_cur(:, i) - r_lag(:, i));
        end

        R_ast = mean(r_lag(:, I), 2) + iA*R_ast*(1/(h*dt));
    end

    return;
end