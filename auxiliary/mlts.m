function [beta, sigma, dres, betaR, sigmaR, dresR] = mlts(x, y, gamma, ns, nc, delta)
% Multivariate Regression
% Input:
%   x: data-matrix (n,p)
%   y: data-matrix (n,q)
%   gamma: proportion of trimming 
%   arguments: 
%     - ns : contains number of subsets; default=5000
%     - nc : number of C-steps; default=10
%     - delta : critical value for Reweighted estimator, deafult=0.01
% Output:
%     beta : matrix (p,q) of MLTS-regression coefficients
%     sigma: matrix (q,q) containing MLTS-residual covariance
%     dres : residual distances (n,1) w.r.t. initial fit
%     betaR : matrix (p,q) of RMLTS-regression coefficients
%     sigmaR: matrix (q,q) containing RMLTS-residual covariance
%     dresR : residual distances (n,1) w.r.t. RMLTS
% Remark: if intercept needed, add a column of ones to the x-matrix
%
% Ref: Agullo,J., Croux, C., and Van Aelst, S. (2008) 
%      The Multivariate Least Trimmed Squares Estimator, 
%      Journal of multivariate analysis
%
% Author: Kristel Joossens

%example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n = 10;
%p = 5;
%q = 2;
%x = [ones(n, 1), randn(n, p-1)]; 
%beta = [2 1 4 2 1; 5 2 1 0.5 2];
%y = x * beta' + randn(n, q);
%gamma = 0.25;
%[beta, sigma, dres, betaR, sigmaR, dresR] = mlts(x, y, gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Added by Michael Pokojovy on 9/15/2024
%% BEGIN
function y = logdet(A)
    E = eig(A);
    y = sum(log(E));
    return;
end
%% END

if nargin < 4
    ns = 500;
end
if nargin < 5
    nc = 10;
end
if nargin < 6
    delta = 0.01;
end

[n, p] = size(x);
q = size(y, 2);
h = floor(n * (1 - gamma)) + 1;
obj0 = 1E10;

for i = 1:ns
    [~, istart] = sort(rand(n, 1));
    istart = istart(1:(p+q));
    xstart = x(istart, :);
    ystart = y(istart, :);
    bstart = (xstart' * xstart) \ (xstart' * ystart);
    sigmastart = (ystart - xstart * bstart)' * (ystart - xstart * bstart) / q;
    
    for j = 1:nc
        res = y - x * bstart;
        tres = res';
        dist2 = sum(pinv(sigmastart) * tres .* tres, 1) / q;
        [~, idist2] = sort(dist2);
        idist2 = idist2(1:h);
        xstart = x(idist2, :);
        ystart = y(idist2, :);
        bstart = (xstart' * xstart) \ (xstart' * ystart);
        sigmastart = (ystart - xstart * bstart)' * (ystart - xstart * bstart) / (h - p);
    end
    
    %obj = det(sigmastart);
    obj = logdet(sigmastart); % replaced det with logdet -- MP 9/15/2024
    if obj < obj0
        result.beta = bstart;
        result.sigma = sigmastart;
        obj0 = obj;
    end
end

cgamma = (1 - gamma) / chi2cdf(chi2inv(1 - gamma, q), q + 2);
result.sigma = cgamma * result.sigma;
res = y - x * result.beta;
tres = res';
result.dres = sqrt(sum(pinv(result.sigma) * tres .* tres, 1)); 

qdelta = sqrt(chi2inv(1 - delta, q));
good = result.dres <= qdelta;
xgood = x(good, :);
ygood = y(good, :);
result.betaR = (xgood' * xgood) \ (xgood' * ygood);
result.sigmaR = (ygood - xgood * result.betaR)' * (ygood - xgood * result.betaR) / (sum(good) - p);
cdelta = (1 - delta) / chi2cdf(qdelta^2, q + 2);
result.sigmaR = cdelta * result.sigmaR;
resR = y - x * result.betaR;
tresR = resR';
result.dresR = sqrt(sum(pinv(result.sigmaR) * tresR .* tresR, 1)); 

beta = result.beta;
sigma = result.sigma;
dres = result.dres;
betaR = result.betaR;
sigmaR = result.sigmaR;
dresR = result.dresR;
end

