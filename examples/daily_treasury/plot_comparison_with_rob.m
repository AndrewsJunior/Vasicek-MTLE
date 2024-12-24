function plot_comparison_with_rob(time_grid_historic, time_grid_backtest, ...
                                  US_rates_historic, EU_rates_historic, ...
                                  US_rates_backtest, EU_rates_backtest, ...
                                  projection_ML, projection_OLS, projection_mlts, projection_mtle, conf_level, OutlierFraction)


    fig = figure(2);

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 26; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);

    time_grid = [time_grid_historic time_grid_backtest(2:end)];
        
    US_rates = [US_rates_historic US_rates_backtest(2:end)];
    EU_rates = [EU_rates_historic EU_rates_backtest(2:end)];
    
    low_rate = 0.50*(min([US_rates EU_rates]) + max([US_rates EU_rates])) - ...
               0.75*(max([US_rates EU_rates]) - min([US_rates EU_rates]));

    high_rate = 0.50*(min([US_rates EU_rates]) + max([US_rates EU_rates])) + ...
                1.00*(max([US_rates EU_rates]) - min([US_rates EU_rates]));


    % Adjust padding values in subplot_tight for better spacing
    padding = [0.12 0.05]; % Increase spacing between subplots

    %% Projection 1
    subplot_tight(2, 2, 3, padding);
    %hold on;
    
    US_pred = plot(time_grid_backtest, mean(squeeze(projection_ML(1, :, :)), 2), ...
                   'Color', [0,      0.4470, 0.7410], 'LineStyle', '--');
    hold on;
    EU_pred = plot(time_grid_backtest, mean(squeeze(projection_ML(2, :, :)), 2), ...
                   'Color', [0.8500, 0.3250, 0.0980], 'LineStyle', '--');
    
    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_ML(1, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_ML(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0,      0.4470, 0.7410], 'FaceAlpha', 0.1);

    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_ML(2, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_ML(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.1);

    plot(time_grid_backtest, quantile(squeeze(projection_ML(1, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');

    plot(time_grid_backtest, quantile(squeeze(projection_ML(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    plot(time_grid_backtest, quantile(squeeze(projection_ML(2, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_ML(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    US = plot(time_grid, US_rates, 'Color', [0,      0.4470, 0.7410]);
    EU = plot(time_grid, EU_rates, 'Color', [0.8500, 0.3250, 0.0980]);
    plot(time_grid_backtest([1 1]), [-10, 10], 'k:');
    
    xlabel("Date", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("Daily short rate", 'interpreter', 'latex', 'FontSize', 8);
    xlim([min(time_grid) max(time_grid)]);
    ylim([low_rate (high_rate + 2)]);

    title('MLE', 'interpreter', 'latex', 'FontSize', 8);
    
    legend([US, EU, US_pred, EU_pred], {'US obs', 'EU obs', 'US mean proj', 'EU mean proj'}, ...
           'Location', 'NorthWest', 'interpreter', 'latex', 'FontSize', 7);

    %% Projection 2
    subplot_tight(2, 2, 4, padding);
    %hold on;

    US_pred = plot(time_grid_backtest, mean(squeeze(projection_ML(1, :, :)), 2), ...
                   'Color', [0,      0.4470, 0.7410], 'LineStyle', '--');
    hold on;
    EU_pred = plot(time_grid_backtest, mean(squeeze(projection_ML(2, :, :)), 2), ...
                   'Color', [0.8500, 0.3250, 0.0980], 'LineStyle', '--');
    
    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_OLS(1, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_OLS(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0,      0.4470, 0.7410], 'FaceAlpha', 0.1);
    
    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_OLS(2, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_OLS(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.1);
    


    plot(time_grid_backtest, quantile(squeeze(projection_OLS(1, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_OLS(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    plot(time_grid_backtest, quantile(squeeze(projection_OLS(2, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_OLS(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    US = plot(time_grid, US_rates, 'Color', [0,      0.4470, 0.7410]);
    EU = plot(time_grid, EU_rates, 'Color', [0.8500, 0.3250, 0.0980]);
    plot(time_grid_backtest([1 1]), [-10, 10], 'k:');
    
    xlabel("Date", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("Daily short rate", 'interpreter', 'latex', 'FontSize', 8);
    xlim([min(time_grid) max(time_grid)]);
    ylim([low_rate (high_rate + 2)]);
    
    title('MLS', 'interpreter', 'latex', 'FontSize', 8);

    legend([US, EU, US_pred, EU_pred], {'US obs', 'EU obs', 'US mean proj', 'EU mean proj'}, ...
           'Location', 'NorthWest', 'interpreter', 'latex', 'FontSize', 7);

    
    %% Projection 3
    subplot_tight(2, 2, 2, padding);
    %hold on;


    US_pred = plot(time_grid_backtest, mean(squeeze(projection_mlts(1, :, :)), 2), ...
                   'Color', [0,      0.4470, 0.7410], 'LineStyle', '--');
    hold on;
    EU_pred = plot(time_grid_backtest, mean(squeeze(projection_mlts(2, :, :)), 2), ...
                   'Color', [0.8500, 0.3250, 0.0980], 'LineStyle', '--');

    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_mlts(1, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_mlts(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0,      0.4470, 0.7410], 'FaceAlpha', 0.1);
    
    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_mlts(2, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_mlts(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.1);
    
    

    plot(time_grid_backtest, quantile(squeeze(projection_mlts(1, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_mlts(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    plot(time_grid_backtest, quantile(squeeze(projection_mlts(2, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_mlts(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');

    
    US = plot(time_grid, US_rates, 'Color', [0,      0.4470, 0.7410]);
    EU = plot(time_grid, EU_rates, 'Color', [0.8500, 0.3250, 0.0980]);
    plot(time_grid_backtest([1 1]), [-10, 10], 'k:');
    
    xlabel("Date", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("Daily short rate", 'interpreter', 'latex', 'FontSize', 8);
    xlim([min(time_grid) max(time_grid)]);
    ylim([low_rate (high_rate + 2)]);

    title(['MLTS ($$\alpha = ' num2str(OutlierFraction) '$$)'], 'interpreter', 'latex', 'FontSize', 8);
    
    legend([US, EU, US_pred, EU_pred], {'US obs', 'EU obs', 'US mean proj', 'EU mean proj'}, ...
           'Location', 'NorthWest', 'interpreter', 'latex', 'FontSize', 7);


        %% Projection 4
    subplot_tight(2, 2, 1, padding);
    %hold on;

    US_pred = plot(time_grid_backtest, mean(squeeze(projection_mtle(1, :, :)), 2), ...
                   'Color', [0,      0.4470, 0.7410], 'LineStyle', '--');
    hold on;
    EU_pred = plot(time_grid_backtest, mean(squeeze(projection_mtle(2, :, :)), 2), ...
                   'Color', [0.8500, 0.3250, 0.0980], 'LineStyle', '--');
    
    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_mtle(1, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_mtle(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0,      0.4470, 0.7410], 'FaceAlpha', 0.1);
    
    patch([time_grid_backtest, flip(time_grid_backtest)], ...
          [quantile(squeeze(projection_mtle(2, :, :)), 0.5*(1.0 - conf_level), 2)', ...
          flip(quantile(squeeze(projection_mtle(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2))'], [0.9, 0.9, 0.9], ...
          'EdgeColor', 'none', 'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.1);
    
    

    plot(time_grid_backtest, quantile(squeeze(projection_mtle(1, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_mtle(1, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    plot(time_grid_backtest, quantile(squeeze(projection_mtle(2, :, :)), 0.5*(1.0 - conf_level), 2), 'k:');
    plot(time_grid_backtest, quantile(squeeze(projection_mtle(2, :, :)), 1.0 - 0.5*(1.0 - conf_level), 2), 'k:');
    
    US = plot(time_grid, US_rates, 'Color', [0,      0.4470, 0.7410]);
    EU = plot(time_grid, EU_rates, 'Color', [0.8500, 0.3250, 0.0980]);
    plot(time_grid_backtest([1 1]), [-10, 10], 'k:');
    
    xlabel("Date", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("Daily short rate", 'interpreter', 'latex', 'FontSize', 8);
    xlim([min(time_grid) max(time_grid)]);
    ylim([low_rate (high_rate + 2)]);

    title(['MTLE ($$\alpha = ' num2str(OutlierFraction) '$$)'], 'interpreter', 'latex', 'FontSize', 8);
    
    legend([US, EU, US_pred, EU_pred], {'US obs', 'EU obs', 'US mean proj', 'EU mean proj'}, ...
           'Location', 'NorthWest', 'interpreter', 'latex', 'FontSize', 7);

    fig.Position = [100, 100, 1200, 600];%800
    print(fig, 'USEUproj1.eps', '-depsc', '-r300');
    hold off;


end
