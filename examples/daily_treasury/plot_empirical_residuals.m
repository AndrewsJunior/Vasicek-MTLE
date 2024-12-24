function plot_empirical_residuals(US_rates, EU_rates, dt, alpha, ...
                                  R_ast_0, A_0, Sigma_0, ...
                                  R_ast_1, A_1, Sigma_1,R_ast_2,A_2, Sigma_2,R_ast_3, A_3, Sigma_3)
    function plot_circle(radius)
        phi_array = linspace(0, 2*pi, 501);
        plot(radius*cos(phi_array(1:end - 1)), radius*sin(phi_array(1:end - 1)), 'k--', 'LineWidth', 2)
    end

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 24; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);

    fig = figure(4);

    n = length(US_rates) - 1;

    % Adjust padding values in subplot_tight for better spacing
    padding = [0.32 0.05]; % Increase spacing between subplots


    %% MLTS
    subplot_tight(1, 4, 2, padding);
    %hold on

    %axis equal;
    %axis([-10 10 -10 10]);

    dr = [diff(US_rates, 1); diff(EU_rates, 1)]' - ...
         (R_ast_2 - [US_rates(1:end - 1); EU_rates(1:end - 1)])'*A_2'*dt;

    res = dr*inv(sqrtm(Sigma_2)*dt);
    plot(res(:, 1), res(:, 2), '*', 'Color', [0.8500 0.3250 0.0980]);
    hold on;

    alpha_cor = 1 - (1 - alpha)^(1/n);
    radius = sqrt(chi2inv(1 - alpha_cor, 2));
    plot_circle(radius);

    I = find(sum(res.*res, 2) >= radius^2);
    plot(res(I, 1), res(I, 2), '*', 'Color', [0.8500 0.3250 0.0980]);%bluedeep[0 0.4470 0.7410]
    axis([-10 10 -10 10]);

    title('MLTS', 'interpreter', 'latex', 'FontSize', 7);

    xlabel("US residuals", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("EU residuals", 'interpreter', 'latex', 'FontSize', 8);





    %% MTLE
    subplot_tight(1, 4, 1, padding);
    %hold on

    %axis equal;
    %axis([-10 10 -10 10]); 

    dr = [diff(US_rates, 1); diff(EU_rates, 1)]' - ...
         (R_ast_3 - [US_rates(1:end - 1); EU_rates(1:end - 1)])'*A_3'*dt;

    res = dr*inv(sqrtm(Sigma_3)*dt);
    plot(res(:, 1), res(:, 2), '*', 'Color', [0.3010 0.7450 0.9330]);%[0.8500 0.3250 0.0980]orange
    hold on;

    alpha_cor = 1 - (1 - alpha)^(1/n);
    radius = sqrt(chi2inv(1 - alpha_cor, 2));
    plot_circle(radius);

    I = find(sum(res.*res, 2) >= radius^2);
    plot(res(I, 1), res(I, 2), '*', 'Color', [0.3010 0.7450 0.9330]);
    axis([-10 10 -10 10]);

    title('MTLE', 'interpreter', 'latex', 'FontSize', 7);

    xlabel("US residuals", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("EU residuals", 'interpreter', 'latex', 'FontSize', 8);



    %% MLE
    subplot_tight(1, 4, 3, padding);
    %hold on

    %axis equal;
    %axis([-10 10 -10 10]);

    dr = [diff(US_rates, 1); diff(EU_rates, 1)]' - ...
         (R_ast_0 - [US_rates(1:end - 1); EU_rates(1:end - 1)])'*A_0'*dt;

    res = dr*inv(sqrtm(Sigma_0)*dt);
    plot(res(:, 1), res(:, 2), '*', 'Color', [0.4940 0.1840 0.5560]);
    hold on;

    alpha_cor = 1 - (1 - alpha)^(1/n);
    radius = sqrt(chi2inv(1 - alpha_cor, 2));
    plot_circle(radius);

    I = find(sum(res.*res, 2) >= radius^2);
    plot(res(I, 1), res(I, 2), '*', 'Color', [0.4940 0.1840 0.5560]);
    axis([-10 10 -10 10]);

    title('MLE', 'interpreter', 'latex', 'FontSize', 7);

    xlabel("US residuals", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("EU residuals", 'interpreter', 'latex', 'FontSize', 8);

    %% OLS
    subplot_tight(1, 4, 4, padding);
    %hold on

    %axis equal;
    %axis([-10 10 -10 10]);

    dr = [diff(US_rates, 1); diff(EU_rates, 1)]' - ...
         (R_ast_1 - [US_rates(1:end - 1); EU_rates(1:end - 1)])'*A_1'*dt;

    res = dr*inv(sqrtm(Sigma_1)*dt);
    plot(res(:, 1), res(:, 2),'*', 'Color', [0 0.4470 0.7410]);
    hold on;

    alpha_cor = 1 - (1 - alpha)^(1/n);
    radius = sqrt(chi2inv(1 - alpha_cor, 2));
    plot_circle(radius);

    I = find(sum(res.*res, 2) >= radius^2);
    plot(res(I, 1), res(I, 2),'*', 'Color', [0 0.4470 0.7410]);
    axis([-10 10 -10 10]);

    title('MLS', 'interpreter', 'latex', 'FontSize', 7);

    xlabel("US residuals", 'interpreter', 'latex', 'FontSize', 8);
    ylabel("EU residuals", 'interpreter', 'latex', 'FontSize', 8);


    fig.Position = [100, 100, 1200, 600];
    print(fig, 'USEUresid1.eps', '-depsc', '-r300');
    hold off;


end