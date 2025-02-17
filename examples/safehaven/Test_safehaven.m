auxiliary = "../../auxiliary";

addpath(auxiliary);

file_name = '???.csv';

opts = detectImportOptions(file_name);

tab = readtable(file_name, opts);
tab{:, 2:7} = log(tab{:, 2:7});
assets = {'XAU', 'XAG', 'BRE',  'WTI', 'CHF', 'JPY'};

numassets = size(tab, 2) - 1;
numrec = size(tab, 1);

incr = 50;
T_array = incr:numrec;

R_ast_hat_0 = zeros(length(T_array), numassets); % OLS
R_ast_hat_1 = zeros(length(T_array), numassets); % MLE
R_ast_hat_2 = zeros(length(T_array), numassets); % MLTS
R_ast_hat_3 = zeros(length(T_array), numassets); % MTLE

dt = 1;
bdpval = 0.25;
rng(1);

parfor i_T = 1:length(T_array)
    R_ast_hat_0(i_T,:) = inverse_map(tab{(i_T):T_array(i_T),2:7}, dt, "ols",  1E-6, 5000);
    R_ast_hat_1(i_T,:) = inverse_map(tab{(i_T):T_array(i_T),2:7}, dt, "ml",   1E-6, 5000);

    R_ast_hat_2(i_T,:) = inverse_map_rob(tab{(i_T):T_array(i_T),2:7}, dt, bdpval, "mlts", 1E-6, 5000);
    R_ast_hat_3(i_T,:) = inverse_map_rob(tab{(i_T):T_array(i_T),2:7}, dt, bdpval, "mtle", 1E-6, 5000);
end

%start timing
tic;

fig = figure(1);
 % Adjust padding values in subplot_tight for better spacing
 padding = [0.09 0.05]; % Increase spacing between subplots

for i_plot = 1:numassets
    subplot_tight(2, 3, i_plot, padding);
    plot(T_array, R_ast_hat_3(:,i_plot), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);%MTLE
    hold on;
    plot(T_array, R_ast_hat_2(:,i_plot), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.','LineWidth', 1.5);%MLTS
    plot(T_array, R_ast_hat_1(:,i_plot), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);%MLE
    plot(T_array, R_ast_hat_0(:,i_plot), 'Color', [0 0.4470 0.7410], 'LineStyle','--','LineWidth', 1.5);%OLS
    
    

    ylim([-3 3]);
    
    title(assets(i_plot), 'interpreter', 'latex', 'FontSize', 14);
    xlabel("Subset index", 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Estimator $\hat{\bf{R}}^{\ast}$', 'interpreter', 'latex', 'FontSize', 18);
    legend('MTLE', 'MLTS', 'MLE', 'MLS',  'Location', 'Best', 'interpreter', 'latex', 'FontSize', 8);
end

fig.Position = [100, 100, 800, 600]; 
figname = 'assets';
print(fig, [char(figname),'.eps'], '-depsc', '-r300');
hold off;

elapsedtime = toc;
fprintf('Elapsed time for simulation: %.2f seconds\n', elapsedtime);

rmpath(auxiliary);

