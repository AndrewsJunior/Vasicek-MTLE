function makeplots_rmse2D(R_all, A_all, Sigma_all, figname_input)
    

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 30; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);

    T_array = 60: 30: 360;

    fig = figure(1);

    hold on;
    subplot_tight(1, 3, 1, [0.25 0.06]);

    ymax = min(25.0, max(R_all,[],"all"))*1.25;
    xlim([min(T_array) max(T_array)]);
    ylim([0 ymax]);

    plot(T_array, R_all(:,1), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);
    hold on;
    plot(T_array, R_all(:,2), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
    plot(T_array, R_all(:,3), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);
    plot(T_array, R_all(:,4), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5);
    xlim([min(T_array) max(T_array)]);
    ylim([0 ymax*1.1]);

    title('Estimator $\hat{\bf{R}}^{\ast}$', 'interpreter', 'latex', 'FontSize', 10);
    xlabel("Calibration time horizon (days)", 'interpreter', 'latex', 'FontSize', 10);
    ylabel("Empirical root-MSE", 'interpreter', 'latex', 'FontSize', 10);
    legend('MTLE', 'MLTS','MLE', 'MLS', 'Location', 'NorthEast', 'interpreter', 'latex', 'FontSize', 8);

    hold on;
    subplot_tight(1, 3, 2, [0.25 0.06]);
    %hold on;
    ymax = min(25.0, max(A_all,[],"all"))*1.25;
    xlim([min(T_array) max(T_array)]);
    ylim([0 ymax]);

    plot(T_array, A_all(:,1), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);
    hold on;
    plot(T_array, A_all(:,2), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
    plot(T_array, A_all(:,3), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);
    plot(T_array, A_all(:,4), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5);
    xlim([min(T_array) max(T_array)]);
    %ylim([0 ymax*1.1]);

    title('Estimator $\hat{\bf{A}}$', 'interpreter', 'latex', 'FontSize', 10);
    xlabel("Calibration time horizon (days)", 'interpreter', 'latex', 'FontSize', 10);
    ylabel("Empirical root-MSE", 'interpreter', 'latex', 'FontSize', 10);
    legend('MTLE', 'MLTS','MLE', 'MLS', 'Location', 'NorthEast', 'interpreter', 'latex', 'FontSize', 8);

    hold on;
    subplot_tight(1, 3, 3, [0.25 0.06]);
    %hold on;
    ymax = max(Sigma_all, [], "all")*1.25;
    xlim([min(T_array) max(T_array)]);
    ylim([0 ymax]);

    plot(T_array, Sigma_all(:,1), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);
    hold on;
    plot(T_array, Sigma_all(:,2), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
    plot(T_array, Sigma_all(:,3), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);
    plot(T_array, Sigma_all(:,4), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5); 


    xlim([min(T_array) max(T_array)]);
    %ylim([0 ymax*1.1]);
    
    title('Estimator $\hat{\bf{\Sigma}}$', 'interpreter', 'latex', 'FontSize', 10);
    xlabel("Calibration time horizon (days)", 'interpreter', 'latex', 'FontSize', 10);
    ylabel("Empirical root-MSE", 'interpreter', 'latex', 'FontSize', 10);
    legend('MTLE', 'MLTS','MLE', 'MLS', 'Location', 'NorthEast', 'interpreter', 'latex', 'FontSize', 8);


    if not(isfolder("SimResultsPlots"))
        mkdir("SimResultsPlots")
    end

    cd ./SimResultsPlots/
    fig.Position = [1, 49, 1920, 955];
    figname = figname_input;
    print(fig, [char(figname),'.eps'], '-depsc', '-r300');

    cd ..\

    close(fig);

end