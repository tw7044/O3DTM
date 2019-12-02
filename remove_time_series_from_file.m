function target_path = remove_time_series_from_file(file_path, target_path)
% REMOVE_TIME_SERIES_FROM_FILE removes time series data from a temperature
% data file and saves this smaller copy for use elsewhere (as time series
% data can be >1gb in size)
%
% REMOVE_TIME_SERIES_FROM_FILE removes time series data from all simulation
% files
% 
% target_path = REMOVE_TIME_SERIES_FROM_FILE(file_path) removes time series
% data from specified file and returns path data is saved to
% 
% target_path = REMOVE_TIME_SERIES_FROM_FILE(file_path, target_path) saves
% data at specified subfolder

if nargin == 0
    source_folder = create_static_path('outputs/ls_temperatures');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            remove_time_series_from_file(file_path);
        end
    end
    
    source_folder = create_static_path('outputs/ls_temperatures/parameters');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            remove_time_series_from_file(file_path, 'parameters');
        end
    end
    
    source_folder = create_static_path('outputs/ls_temperatures/diviner_fit');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            remove_time_series_from_file(file_path, 'diviner_fit');
        end
    end
    fprintf('ALL DONE\n')
    return
end
target_folder = create_static_path('outputs/ls_temperatures_extremes');

if nargin > 1
    target_folder = sprintf('%s/%s', target_folder, target_path);
end

file_name = strsplit(file_path, '/');
file_name = file_name{end};

fprintf('%s ', file_name)

data = load(file_path);
data = data.data;
data.time_series_extremes = struct;
data.time_series_extremes.description = 'Values calculated from time series data (full time series data not included in this file to save space)';
data.time_series_extremes.Tmax_3dmat = max(data.time_series.T_time_series_4dmat(:,:,:,:), [], 4);
data.time_series_extremes.Tmin_3dmat = min(data.time_series.T_time_series_4dmat(:,:,:,:), [], 4);
data.time_series_extremes.Tmean_3dmat = mean(data.time_series.T_time_series_4dmat(:,:,:,:), 4);
data = rmfield(data, 'time_series');

target_name = sprintf('%s_extremes.mat', strrep(file_name, '.mat', ''));
target_path = sprintf('%s/%s', target_folder, target_name);
save(target_path, 'data')

fprintf('DONE\n')
end