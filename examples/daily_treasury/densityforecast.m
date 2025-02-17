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
%figure(1);
%plot_short_rates(time_grid_historic, time_grid_backtest, ...
                 %US_rates_historic, EU_rates_historic, ...
                 %US_rates_backtest, EU_rates_backtest);

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

nrep = 5000;%10000;
n = length(time_grid_backtest);

pred_var  = forward_map(R_ast_var,  A_var,  Sigma_var,  R0, dt, n, nrep);%MLE
pred_ols  = forward_map(R_ast_ols,  A_ols,  Sigma_ols,  R0, dt, n, nrep);%OLS
pred_mlts = forward_map(R_ast_mlts, A_mlts, Sigma_mlts, R0, dt, n, nrep);%MLTS
pred_mtle = forward_map(R_ast_mtle, A_mtle, Sigma_mtle, R0, dt, n, nrep);%MTLE

fig = figure(1);
% Adjust padding values in subplot_tight for better spacing
padding = [0.32 0.05]; % Increase spacing between subplots

%contour
xgrid = linspace(5.1, 6);
ygrid = linspace(3, 4.4);
[x1, x2] = meshgrid(xgrid, ygrid);
xi = [x1(:) x2(:)];

%MTLE
subplot_tight(1, 4, 1, padding);

MTLEUS = squeeze(pred_mtle(1, :, :));
MTLEUSPREDS = MTLEUS(end, :);

MTLEEU = squeeze(pred_mtle(2, :, :));
MTLEEUPREDS = MTLEEU(end, :);

[f4, ep4] = ksdensity([MTLEUSPREDS' MTLEEUPREDS'], xi);

X = reshape(ep4(:,1), length(xgrid), length(ygrid));
Y = reshape(ep4(:,2), length(xgrid), length(ygrid));
Z = reshape(f4, length(xgrid), length(ygrid));
%figure
contour(X,Y,Z,10)

title('MTLE', 'interpreter', 'latex', 'FontSize', 7);

%MLTS
subplot_tight(1, 4, 2, padding);

MLTSUS = squeeze(pred_mlts(1, :, :));
MLTSUSPREDS = MLTSUS(end, :);

MLTSEU = squeeze(pred_mlts(2, :, :));
MLTSEUPREDS = MLTSEU(end, :);

[f3, ep3] = ksdensity([MLTSUSPREDS' MLTSEUPREDS'], xi);

X = reshape(ep3(:,1), length(xgrid), length(ygrid));
Y = reshape(ep3(:,2), length(xgrid), length(ygrid));
Z = reshape(f3, length(xgrid), length(ygrid));
%figure
contour(X,Y,Z,10)

title('MLTS', 'interpreter', 'latex', 'FontSize', 7);


%MLE
subplot_tight(1, 4, 3, padding);

MLEUS = squeeze(pred_var(1, :, :));
MLEUSPREDS = MLEUS(end, :);

MLEEU = squeeze(pred_var(2, :, :));
MLEEUPREDS = MLEEU(end, :);

[f1,ep1] = ksdensity([MLEUSPREDS' MLEEUPREDS'], xi);

X = reshape(ep1(:,1), length(xgrid), length(ygrid));
Y = reshape(ep1(:,2), length(xgrid), length(ygrid));
Z = reshape(f1, length(xgrid), length(ygrid));
contour(X,Y,Z,10)

title('MLE', 'interpreter', 'latex', 'FontSize', 7);

%OLS
subplot_tight(1, 4, 4, padding);

OLSUS = squeeze(pred_ols(1, :, :));
OLSUSPREDS = OLSUS(end, :);

OLSEU = squeeze(pred_ols(2, :, :));
OLSEUPREDS = OLSEU(end, :);

[f2, ep2] = ksdensity([OLSUSPREDS' OLSEUPREDS'], xi);

X = reshape(ep2(:,1), length(xgrid), length(ygrid));
Y = reshape(ep2(:,2), length(xgrid), length(ygrid));
Z = reshape(f2, length(xgrid), length(ygrid));
%figure
contour(X,Y,Z,10)

title('MLS', 'interpreter', 'latex', 'FontSize', 7);

fig.Position = [100, 100, 1200, 600];
print(fig, 'USEUdensity.eps', '-depsc', '-r300');
%hold off;

rmpath(auxiliary);
