function plot_short_rates(time_grid_historic, time_grid_backtest, ...
                          US_rates_historic, EU_rates_historic, ...
                          US_rates_backtest, EU_rates_backtest)

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
                0.75*(max([US_rates EU_rates]) - min([US_rates EU_rates]));

    % Rate vs Date
    subplot_tight(1, 2, 1, [0.08 0.06]);
    hold on

    title('3 Mo Yields vs Date', 'interpreter', 'latex', 'FontSize', 18);

    plot(time_grid, US_rates);
    plot(time_grid, EU_rates);
    plot(time_grid_backtest([1 1]), [-10, 10], 'k:');

    xlabel("Date", 'interpreter', 'latex', 'FontSize', 18);
    ylabel("Daily rate", 'interpreter', 'latex', 'FontSize', 18);
    xlim([min(time_grid) max(time_grid)]);
    ylim([low_rate high_rate]);

    legend('US', 'EU', 'Location', 'NorthWest', 'interpreter', 'latex', 'FontSize', 18);

    % EU rate vs US rate
    subplot_tight(1, 2, 2, [0.08 0.06]);
    hold on

    title('US vs EU 3 Mo Yield Differences', 'interpreter', 'latex', 'FontSize', 18);

    plot(diff(US_rates), diff(EU_rates), '*');
    plot(mean(diff(US_rates)), mean(diff(EU_rates)), 'r.', 'MarkerSize', 25);
    xlabel("Daily US rate increment", 'interpreter', 'latex', 'FontSize', 18);
    ylabel("Daily EU rate increment", 'interpreter', 'latex', 'FontSize', 18);

    hold off
end