function plot_error_backtest(time_grid_backtest, US_rates_backtest, EU_rates_backtest, ...
                             pred_0, pred_1, pred_2, pred_3)

    %0: MLE
    %1: OLS
    %2: MLTS
    %3: MTLE

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 26; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);

    fig = figure(3);

    [~, n, nrep] = size(pred_0);

    MSE0 = zeros(n, 1);
    MSE1 = zeros(n, 1);
    MSE2 = zeros(n, 1);
    MSE3 = zeros(n, 1);

    MAPE0 = zeros(n, 1);
    MAPE1 = zeros(n, 1);
    MAPE2 = zeros(n, 1);
    MAPE3 = zeros(n, 1);

    for i = 1:n
    for rep = 1:nrep
        x = [US_rates_backtest(i); EU_rates_backtest(i)];

        MSE0(i) = MSE0(i) + sum((pred_0(:, i, rep) - x).^2)/nrep;%MLE
        MSE1(i) = MSE1(i) + sum((pred_1(:, i, rep) - x).^2)/nrep;%OLS
        MSE2(i) = MSE2(i) + sum((pred_2(:, i, rep) - x).^2)/nrep;%MLTS
        MSE3(i) = MSE3(i) + sum((pred_3(:, i, rep) - x).^2)/nrep;%MTLE
            
        MAPE0(i) = MAPE0(i) + 100.0*sum(abs((pred_0(:, i, rep) - x)./x))/nrep;%MLE
        MAPE1(i) = MAPE1(i) + 100.0*sum(abs((pred_1(:, i, rep) - x)./x))/nrep;%OLS
        MAPE2(i) = MAPE2(i) + 100.0*sum(abs((pred_2(:, i, rep) - x)./x))/nrep;%MLTS
        MAPE3(i) = MAPE3(i) + 100.0*sum(abs((pred_3(:, i, rep) - x)./x))/nrep;%MTLE
    end
    end

    % root-MSE vs Date
    
    padding = [0.22 0.07];
    subplot_tight(1, 2, 1, padding);

    % Using RGB triplets for specific colors

    plot(time_grid_backtest, sqrt(MSE3), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);  % (?)Dark yellow
    hold on;
    plot(time_grid_backtest, sqrt(MSE2), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);  % Dark orange
    plot(time_grid_backtest, sqrt(MSE0), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5);  % Dark purple
    plot(time_grid_backtest, sqrt(MSE1), 'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5);  % Dark blue

    xlabel("Date", 'interpreter', 'latex', 'FontSize', 18);
    ylabel("Root-MSE", 'interpreter', 'latex', 'FontSize', 18);
    title('Empirical Root-MSE vs Date', 'interpreter', 'latex', 'FontSize', 18);

    legend({'MTLE','MLTS','MLE', 'MLS'}, ...
           'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);

    % MPE vs Date
    subplot_tight(1, 2, 2, padding);
    %hold on
    
    plot(time_grid_backtest, MAPE3, 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5);  
    hold on;
    plot(time_grid_backtest, MAPE2, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', 1.5);
    plot(time_grid_backtest, MAPE0, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', ':', 'LineWidth', 1.5); 
    plot(time_grid_backtest, MAPE1,  'Color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth', 1.5);
    
    

    xlabel("Date", 'interpreter', 'latex', 'FontSize', 18);
    ylabel("MAPE (\%)", 'interpreter', 'latex', 'FontSize', 18);
    title('Empirical MAPE vs Date', 'interpreter', 'latex', 'FontSize', 18);

    legend({'MTLE','MLTS','MLE', 'MLS'}, ...
           'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);

    fig.Position = [100, 100, 1200, 600];
    print(fig, 'USEUbacktest1.eps', '-depsc', '-r300');
    hold off;
end