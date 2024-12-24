%Source: https://github.com/mpokojovy/mvvasicek
function sol = forward_map(R_ast, A, Sigma, R0, dt, n, nrep)
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
% Copyright (c) 2023: Michael Pokojovy, Ebenezer Nkum and
%                     Thomas M. Fullerton, Jr.

    p  = length(R_ast);

    a  = A*R_ast*dt;
    b  = eye(p) - A*dt;
    c  = sqrtm(Sigma)*sqrt(dt);
    
    sol = zeros(p, n, nrep);
    
    for rep = 1:nrep
        sol(:, 1 , rep) = R0;
           
        for i = 2:n
            sol(:, i, rep) = a + b*sol(:, i-1, rep) + c*randn(p, 1);
        end 
    end 
end