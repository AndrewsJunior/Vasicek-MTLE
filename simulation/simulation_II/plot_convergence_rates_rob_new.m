function plot_convergence_rates_rob_new(R_ast_errors_0, R_ast_errors_1, R_ast_errors_2, R_ast_errors_3, A_errors_0, A_errors_1, A_errors_2, ...
    A_errors_3, Sigma_errors_0, Sigma_errors_1, Sigma_errors_2, Sigma_errors_3, figname_input)

set(gcf, 'PaperUnits', 'centimeters');
xSize = 30; ySize = 12;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf,'Position', [0 0 xSize*50 ySize*50]);

T_array = 60: 30: 360;


fig = figure(1);
hold on;
subplot_tight(1, 3, 1, [0.25 0.06]);

ymax = min(25.0, max([sqrt(mean(R_ast_errors_1, 2)); sqrt(mean(R_ast_errors_0, 2)); sqrt(mean(R_ast_errors_2, 2)); sqrt(mean(R_ast_errors_3, 2))]))*1.25;
xlim([min(T_array) max(T_array)]);
ylim([0 ymax]);

plot(T_array, sqrt(mean(R_ast_errors_3, 2)), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(T_array, sqrt(mean(R_ast_errors_2, 2)), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
plot(T_array, sqrt(mean(R_ast_errors_1, 2)), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);
plot(T_array, sqrt(mean(R_ast_errors_0, 2)), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5);
xlim([min(T_array) max(T_array)]);
ylim([0 ymax*1.1]);

title('Estimator $\hat{\bf{R}}^{\ast}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel("Calibration time horizon (days)", 'interpreter', 'latex', 'FontSize', 10);
ylabel("Empirical root-MSE", 'interpreter', 'latex', 'FontSize', 10);
legend('MTLE', 'MLTS','MLE', 'MLS', 'Location', 'NorthEast', 'interpreter', 'latex', 'FontSize', 8);

hold on;
subplot_tight(1, 3, 2, [0.25 0.06]);
    %hold on;
ymax = min(25.0, max([sqrt(mean(A_errors_1, 2)); sqrt(mean(A_errors_0, 2)) ; sqrt(mean(A_errors_2, 2)) ; sqrt(mean(A_errors_3, 2))]))*1.25;
xlim([min(T_array) max(T_array)]);
ylim([0 ymax]);

plot(T_array, sqrt(mean(A_errors_3, 2)), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(T_array, sqrt(mean(A_errors_2, 2)), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
plot(T_array, sqrt(mean(A_errors_1, 2)), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);
plot(T_array, sqrt(mean(A_errors_0, 2)), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5);
xlim([min(T_array) max(T_array)]);
%ylim([0 ymax*1.1]);

title('Estimator $\hat{\bf{A}}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel("Calibration time horizon (days)", 'interpreter', 'latex', 'FontSize', 10);
ylabel("Empirical root-MSE", 'interpreter', 'latex', 'FontSize', 10);
legend('MTLE', 'MLTS','MLE', 'MLS', 'Location', 'NorthEast', 'interpreter', 'latex', 'FontSize', 8);

hold on;
subplot_tight(1, 3, 3, [0.25 0.06]);
    %hold on;
ymax = max([sqrt(mean(Sigma_errors_1, 2)); sqrt(mean(Sigma_errors_0, 2)); sqrt(mean(Sigma_errors_2, 2)); sqrt(mean(Sigma_errors_3, 2))])*1.25;
xlim([min(T_array) max(T_array)]);
ylim([0 ymax]);

plot(T_array, sqrt(mean(Sigma_errors_3, 2)), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(T_array, sqrt(mean(Sigma_errors_2, 2)), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
plot(T_array, sqrt(mean(Sigma_errors_1, 2)), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);
plot(T_array, sqrt(mean(Sigma_errors_0, 2)), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5); 


xlim([min(T_array) max(T_array)]);
%ylim([0 ymax*1.1]);

title('Estimator $\hat{\bf{\Sigma}}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel("Calibration time horizon (days)", 'interpreter', 'latex', 'FontSize', 10);
ylabel("Empirical root-MSE", 'interpreter', 'latex', 'FontSize', 10);
legend('MTLE', 'MLTS','MLE', 'MLS', 'Location', 'NorthEast', 'interpreter', 'latex', 'FontSize', 8);
%NorthEast

if not(isfolder("SimMSEResultsPlots"))
    mkdir("SimMSEResultsPlots")
end

cd ./SimMSEResultsPlots/
fig.Position = [1, 49, 1920, 955];
figname = figname_input;
print(fig, [char(figname),'.eps'], '-depsc', '-r300');

cd ..\

close(fig);


end