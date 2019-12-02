function save_ls_stability_plot(ls_name)
% SAVE_LS_STABILITY_PLOT saves stability plot from simulated temperatures
%
% SAVE_LS_STABILITY_PLOT saves stability plots for all data files
%
% SAVE_LS_STABILITY_PLOT(ls_name) saves stability plots specified landing
% site or data file path

width = 12;
height = 9;

if nargin == 0
    source_folder = create_static_path('outputs/ls_temperatures_extremes');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            load(file_path, 'data')
            if max(data.metadata.z_arr) > 5 || data.metadata.created < datetime(2018,9,8)
                fprintf('%s *** OLD ***\n', listing(idx).name)
                continue
            end
            fprintf('%s ', listing(idx).name);
            plot_ls_stability(file_path, 'multi')
            
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 width height];
            fig.Units = 'inches';
            fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
            drawnow
            
            image_name = create_static_path(sprintf('plots/ls_stability/%s', strrep(listing(idx).name, '_extremes.mat', '_stability.png')));
            print(image_name, '-dpng', '-r300')
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    
    source_folder = create_static_path('outputs/ls_temperatures_extremes/parameters');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            load(file_path, 'data')
            if max(data.metadata.z_arr) > 5 || data.metadata.created < datetime(2018,9,8)
                fprintf('%s *** OLD ***\n', listing(idx).name)
                continue
            end
            fprintf('%s ', listing(idx).name);
            plot_ls_stability(file_path, 'multi')
            
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 width height];
            fig.Units = 'inches';
            fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
            drawnow
            
            image_name = create_static_path(sprintf('plots/ls_stability/parameters/%s', strrep(listing(idx).name, '_extremes.mat', '_stability.png')));
            print(image_name, '-dpng', '-r300')
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')

    source_folder = create_static_path('outputs/ls_temperatures_extremes/diviner_fit');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            load(file_path, 'data')
            if max(data.metadata.z_arr) > 5 || data.metadata.created < datetime(2018,9,8)
                fprintf('%s *** OLD ***\n', listing(idx).name)
                continue
            end
            fprintf('%s ', listing(idx).name);
            plot_ls_stability(file_path, 'multi')
            
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 width height];
            fig.Units = 'inches';
            fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
            drawnow
            
            image_name = create_static_path(sprintf('plots/ls_stability/diviner_fit/%s', strrep(listing(idx).name, '_extremes.mat', '_stability.png')));
            print(image_name, '-dpng', '-r300')
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    return
end

plot_ls_stability(ls_name, 'multi')
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
image_name = create_static_path(sprintf('plots/ls_stability/%s', strrep(ls_name, '.mat', '_stability.png')));
print(image_name, '-dpng', '-r300')
end