auxiliary = "../../auxiliary";

addpath(auxiliary);

time_grid_historic = datetime("1/1/2023", 'InputFormat', 'MM/dd/uuuu'):datetime("12/31/2023", 'InputFormat', 'MM/dd/uuuu');
time_grid_backtest = datetime("1/1/2024", 'InputFormat', 'MM/dd/uuuu'):datetime("3/31/2024", 'InputFormat', 'MM/dd/uuuu');


run load_data.m

US_rates_historic = interp1(US_dates, US_rates, time_grid_historic, 'linear', 'extrap');
EU_rates_historic = interp1(EU_dates, EU_rates, time_grid_historic, 'linear', 'extrap');

US_rates_backtest = interp1(US_dates, US_rates, time_grid_backtest, 'linear', 'extrap');
EU_rates_backtest = interp1(EU_dates, EU_rates, time_grid_backtest, 'linear', 'extrap');

% Plot data
figure(1);
plot_short_rates(time_grid_historic, time_grid_backtest, ...
                 US_rates_historic, EU_rates_historic, ...
                 US_rates_backtest, EU_rates_backtest);

%% Model calibration
dt = 1.0;

R_hist = [US_rates_historic; EU_rates_historic]';
R_bm   = [US_rates_backtest; EU_rates_backtest]';

R0 = R_bm(1, :);

%For the rob method
OutlierFraction = 0.2;
tol = 1e-8;
MaxIter = 1000;

[R_ast_var,  A_var,  Sigma_var]  = inverse_map(R_hist, dt, "ml", 1E-8, 1E6); %MLE
[R_ast_ols,  A_ols,  Sigma_ols]  = inverse_map(R_hist, dt, "ols", 1E-8, 1E6); %OLS
[R_ast_mlts, A_mlts, Sigma_mlts] = inverse_map_rob(R_hist, dt, OutlierFraction, 'mlts', tol, MaxIter);%MLTS
[R_ast_mtle, A_mtle, Sigma_mtle] = inverse_map_rob(R_hist, dt, OutlierFraction, 'mtle', tol, MaxIter);%MTTE


%% Rate projection
conf_level = 0.90;

nrep = 10000;
n = length(time_grid_backtest);

pred_var  = forward_map(R_ast_var,  A_var,  Sigma_var,  R0, dt, n, nrep);%MLE
pred_ols  = forward_map(R_ast_ols,  A_ols,  Sigma_ols,  R0, dt, n, nrep);%OLS
pred_mlts = forward_map(R_ast_mlts, A_mlts, Sigma_mlts, R0, dt, n, nrep);%MLTS
pred_mtle = forward_map(R_ast_mtle, A_mtle, Sigma_mtle, R0, dt, n, nrep);%MTLE

%figure(2);
plot_comparison_with_rob(time_grid_historic, time_grid_backtest, ...
                US_rates_historic, EU_rates_historic, ...
                US_rates_backtest, EU_rates_backtest, ...
                pred_var, pred_ols, pred_mlts, pred_mtle, conf_level, OutlierFraction);


%% Backtesting - Error comparison
%figure(3);
plot_error_backtest(time_grid_backtest, US_rates_backtest, EU_rates_backtest, ...
                    pred_var, pred_ols, pred_mlts, pred_mtle);

%% Empirical residuals
%figure(4);
plot_empirical_residuals(US_rates_historic, EU_rates_historic, dt, 0.05, ...
                         R_ast_var, A_var, Sigma_var, ...
                         R_ast_ols, A_ols, Sigma_ols, ...
                         R_ast_mlts, A_mlts, Sigma_mlts, ...
                         R_ast_mtle, A_mtle, Sigma_mtle);

%% Markowitz optimal portfolio

return_min = 0.011;

disp("MLE");
[portfolio_weights, mu_return, Sigma_return, Sharpe_ratio, RoR, RoR_backtest] = ...
    Markowitz_portfolio(pred_var, [US_rates_backtest; EU_rates_backtest], US_EU_xchange_rate, return_min, dt);
fprintf('RoR = %s\n', sprintf("%5.4f", 100*RoR));

disp("MLS");
[portfolio_weights, mu_return, Sigma_return, Sharpe_ratio, RoR, RoR_backtest] = ...
    Markowitz_portfolio(pred_ols, [US_rates_backtest; EU_rates_backtest], US_EU_xchange_rate, return_min, dt);
fprintf('RoR = %s\n', sprintf("%5.4f", 100*RoR));

disp("MLTS");
[portfolio_weights, mu_return, Sigma_return, Sharpe_ratio, RoR, RoR_backtest] = ...
    Markowitz_portfolio(pred_mlts, [US_rates_backtest; EU_rates_backtest], US_EU_xchange_rate, return_min, dt);
fprintf('RoR = %s\n', sprintf("%5.4f", 100*RoR));

disp("MLTE");
[portfolio_weights, mu_return, Sigma_return, Sharpe_ratio, RoR, RoR_backtest] = ...
    Markowitz_portfolio(pred_mtle, [US_rates_backtest; EU_rates_backtest], US_EU_xchange_rate, return_min, dt);
fprintf('RoR = %s\n', sprintf("%5.4f", 100*RoR));

rmpath(auxiliary);