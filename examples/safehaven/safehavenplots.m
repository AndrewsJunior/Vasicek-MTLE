file_name = '???.csv';

opts = detectImportOptions(file_name);

tab = readtable(file_name, opts);
tab{:, 2:7} = log(tab{:, 2:7});
assets = {'XAU', 'XAG', 'BRE',  'WTI', 'CHF', 'JPY'};

%subplot_tight
%
% Adjust padding values in subplot_tight for better spacing
padding = [0.10 0.050]; % Increase spacing between subplots

%% Projection 1
fig = figure(1);
hold on;
subplot_tight(2, 3, 1, padding);
plot(tab{:,1}, tab{:,2}, 'k-')
title('XAU', 'interpreter', 'latex', 'FontSize', 8);
subplot_tight(2, 3, 2, padding );
plot(tab{:,1}, tab{:,3}, 'k-')
title('XAG', 'interpreter', 'latex', 'FontSize', 8);
subplot_tight(2, 3, 3, padding);
plot(tab{:,1}, tab{:,4}, 'k-')
title('BRE', 'interpreter', 'latex', 'FontSize', 8);


subplot_tight(2, 3, 4, padding);
plot(tab{:,1}, tab{:,5}, 'k-')
title('WTI', 'interpreter', 'latex', 'FontSize', 8);
subplot_tight(2, 3, 5, padding );
plot(tab{:,1}, tab{:,6}, 'k-')
title('CHF', 'interpreter', 'latex', 'FontSize', 8);
subplot_tight(2, 3, 6, padding);
plot(tab{:,1}, tab{:,7}, 'k-')
title('JPY', 'interpreter', 'latex', 'FontSize', 8);


fig.Position = [100, 100, 800, 600];
print(fig, 'safehavenplots.eps', '-depsc', '-r300');