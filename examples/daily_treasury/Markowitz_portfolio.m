function [portfolio_weights, mu_return, Sigma_return, Sharpe_ratio, RoR, RoR_backtest] = ...
          Markowitz_portfolio(pred, backtest, US_EU_xchange_rate, return_min, dt)
    p    = size(pred, 1);
    nrep = size(pred, 3);

    pred     = pred/100.0; % convert from percent to ratios
    backtest = backtest/100.0;
    dt       = dt/365; % all days, not just trading

    RoR = zeros(p, nrep);

    for rep = 1:nrep
        ror = log(1.0 + dt*squeeze(pred(:, 2:end, rep)));
        RoR(:, rep) = exp(sum(ror, 2)) - 1.0;

        %RoR(:, rep) = exp(dt*sum(squeeze(pred(:, 2:end, rep)), 2)) - 1.0; % use exp() for better numerical stability
    end

    mu_return    = mean(RoR, 2);
    Sigma_return = cov(RoR');

    % Markowitz-style portfolio
    % https://docs.mosek.com/portfolio-cookbook/markowitz.html

    c = @(w) w'*Sigma_return*w;
    A = mu_return';

    opts = optimset('Display', 'off');
    portfolio_weights = fmincon(c, zeros(p, 1), -A, -return_min, ones(1, p), 1.0, zeros(p, 1), ones(p, 1), [], opts);

    Sharpe_ratio = (portfolio_weights'*mu_return - return_min)/ ...
                    sqrt(portfolio_weights'*Sigma_return*portfolio_weights);

    RoR = portfolio_weights'*mu_return; % assuming constant exchange rate

    if (~isempty(US_EU_xchange_rate))
        wealth = 1.0;

        for i = 2:length(US_EU_xchange_rate)
            conv = [1.0; 1.0/US_EU_xchange_rate(i)]./[1.0; 1.0/US_EU_xchange_rate(i-1)];
    
            wealth = wealth*(1.0 + dt*portfolio_weights'*(backtest(:, i).*conv));
        end
    
        RoR_backtest = wealth - 1.0;
    else
        RoR_backtest = [];
    end
end