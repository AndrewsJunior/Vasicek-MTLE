function sol = forward_map_out(R_ast, A, Sigma, R0, dt, n, nrep, out_pcent, ncp)
% forward_map simulates random path of a discrete multivariate Vasicek
%             process
%   sol = forward_map(R_ast, A, Sigma, R0, dt, n, nrep)
%
% Inputs: 
%   R_ast : long-term mean (p x 1)
%   A     : reversion speed matrix (p x p)
%   Sigma : volatility matrix (p x p)
%   R0    : vector of initial rates
%   dt    : time step
%   n     : number of time steps
%   nrep  : number of replications
%
% Outputs:
%     sol :  p x n x nrep tensor of nrep simulated p-variate paths
%
% Copyright (c) 2024: Andrews T. Anum, Michael Pokojovy, Ebenezer Nkum and
%                     Thomas M. Fullerton, Jr.

    p  = length(R_ast);

    a  = A*R_ast*dt;
    b  = eye(p) - A*dt;
    c  = sqrtm(Sigma)*sqrt(dt);
    
    sol = zeros(p, n, nrep);
    numOnes = round(out_pcent*n); %num of outliers
    numZeros = (n) - numOnes;
    randomVector = [zeros(1, numZeros), ones(1, numOnes) ];


    for rep = 1:nrep
        sol(:, 1 , rep) = R0;
        for i = 2:n
            if randomVector(i) == 1
                delta = sqrt(ncp/p)*ones(1,p);
                shift = delta + randn(p, 1)'*eye(p);
                eps = shift';
            else
                eps = randn(p, 1);
            end
            sol(:, i, rep) = a + b*sol(:, i-1, rep) + c*eps;
        end 
    end
    
end