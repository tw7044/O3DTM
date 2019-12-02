function save_ls_diviner_comparison_plot(ls_name)
% SAVE_LS_DIVINER_COMPARISON_PLOT saves landing site comparison summary
% plot for model and diviner temperatures
%
% SAVE_LS_DIVINER_COMPARISON_PLOT saves plots for all data files
%
% SAVE_LS_DIVINER_COMPARISON_PLOT(ls_name) saves plots for specified
% landing site or data file location

width = 15;
height = 9;

if nargin == 0
    source_folder = create_static_path('outputs/ls_temperatures');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    
    clf
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 width height];
    fig.Units = 'inches';
    fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
    
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            if ~strcmp(listing(idx).name(1:2), 'ls')
                fprintf('%s * NOT LS *\n', listing(idx).name)
                continue
            end
            load(strrep(strrep(file_path, 'ls_temperatures', 'ls_temperatures_extremes'), '.mat', '_extremes.mat'), 'data')
            if max(data.metadata.z_arr) > 5
                fprintf('%s *** OLD ***\n', listing(idx).name)
                continue
            end
            fprintf('%s ', listing(idx).name);
            rehash
            save_ls_diviner_comparison_plot(file_path)
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    
    source_folder = create_static_path('outputs/ls_temperatures/parameters');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            if ~strcmp(listing(idx).name(1:2), 'ls')
                fprintf('%s * NOT LS *\n', listing(idx).name)
                continue
            end
            load(strrep(strrep(file_path, 'ls_temperatures', 'ls_temperatures_extremes'), '.mat', '_extremes.mat'), 'data')
            if max(data.metadata.z_arr) > 5
                fprintf('%s *** OLD ***\n', listing(idx).name)
                continue
            end
            fprintf('%s ', listing(idx).name);
            rehash
            save_ls_diviner_comparison_plot(file_path)
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    
    source_folder = create_static_path('outputs/ls_temperatures/diviner_fit');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            if ~strcmp(listing(idx).name(1:2), 'ls')
                fprintf('%s * NOT LS *\n', listing(idx).name)
                continue
            end
            fprintf('%s ', listing(idx).name);
            rehash
            save_ls_diviner_comparison_plot(file_path)
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    return
end

plot_ls_diviner_comparison(ls_name);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 width height];
fig.Units = 'inches';
fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
drawnow

parameters_bool = contains(ls_name, 'parameters');
fit_bool = contains(ls_name, '_diviner_fit');
ls_name = strrep(ls_name, '\', '/');
ls_name = strsplit(ls_name, '/');
ls_name = ls_name{end};
ls_name = strrep(ls_name, '_extremes', '');
ls_name = strrep(ls_name, '_.mat', '.mat');
if parameters_bool
    ls_name = sprintf('parameters/%s', ls_name);
end
if fit_bool
    ls_name = sprintf('diviner_fit/%s', ls_name);
end
image_name = create_static_path(sprintf('plots/ls_diviner_comparison/%s', strrep(ls_name, '.mat', '_diviner_comparison.png')));
print(image_name, '-dpng', '-r300')
end