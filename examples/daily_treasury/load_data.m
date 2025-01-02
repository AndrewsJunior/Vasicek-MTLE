%% Time grids
wrn = warning();
warning off

US_dates = [];
US_rates = [];

EU_dates = [];
EU_rates = [];

for year = 2019:2024
    % US rates
    file_name = ['data/US-daily-rates-', num2str(year), '.csv'];
    opts = detectImportOptions(file_name);
    opts = setvartype(opts, 'Date', 'datetime');
    opts = setvaropts(opts, 'Date', 'InputFormat', 'MM/dd/uuuu');
    tab = readtable(file_name, opts);

    US_dates = [US_dates; tab.('Date')];
    US_rates = [US_rates; tab.('x3M')];

    % EU rates
    file_name = ['data/EU-daily-rates-', num2str(year), '.csv'];
    opts = detectImportOptions(file_name);
    opts = setvartype(opts, 'Date', 'datetime');
    opts = setvaropts(opts, 'Date', 'InputFormat', 'MM/dd/uuuu');
    tab = readtable(file_name, opts);

    EU_dates = [EU_dates; tab.('Date')];
    EU_rates = [EU_rates; tab.('x3M')];
end


%% clean up temporary variables
warning(wrn);
clear("file_name", "opts", "tab", "wrn", "year");