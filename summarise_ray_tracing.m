function summarise_ray_tracing(ls_name, use_seasons)
% SUMMARISE_RAY_TRACING producs statistics about ray tracing to save having
% to load large data file
%
% SUMMARISE_RAY_TRACING summarises ray tracing for all data files
%
% SUMMARISE_RAY_TRACING(ls_name) defines landing site name or data file
% path to summarise
%
% SUMMARISE_RAY_TRACING(ls_name, use_seasons) controls if seasonal ray
% tracing is summarised

if nargin == 0
    for seasons_str = {'no_seasons', 'seasons'}
        source_folder = create_static_path(sprintf('crater_environments/ray_tracing/%s', seasons_str{1}));
        fprintf('Processing all files in %s...\n', source_folder)
        listing = dir(source_folder);
        for idx = 1:numel(listing)
            if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
                file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
                summarise_ray_tracing(file_path);
            end
        end
        fprintf('FOLDER DONE\n')
    end
    fprintf('\n ALL DONE\n')
    return
end
if nargin == 1 && ~(contains(ls_name, '/') || contains(ls_name, '\'))
    summarise_ray_tracing(ls_name, false);
    summarise_ray_tracing(ls_name, true);
    return
end

if contains(ls_name, '/') || contains(ls_name, '\')
    ls_name = strrep(ls_name, '\', '/');
    source_path = ls_name;
else
    if use_seasons
        file_season_str = 'seasons';
    else
        file_season_str = 'no_seasons';
    end
    source_path = create_static_path(sprintf('crater_environments/ray_tracing/%s/%s_ray_tracing.mat', file_season_str, ls_name));
end

file_name = strsplit(source_path, '/');
file_name = sprintf('%s/%s', file_name{end-1:end});
target_path = create_static_path(sprintf('crater_environments/ray_tracing/summary/%s_summary.mat', file_name(1:end-4)));
fprintf('%s... ', file_name);

load(source_path, 'data');
if data.created < datetime(2018,9,8)
    fprintf('*** OLD ***\n')
    return
end
illuminated_time_matrix = 1-sum(isnan(data.theta_3dmat),3)/size(data.theta_3dmat,3);
theta_max_matrix = nanmax(data.theta_3dmat, [], 3);
theta_min_matrix = nanmin(data.theta_3dmat, [], 3);
theta_mean_matrix = nanmean(data.theta_3dmat, 3);

data = rmfield(data, 'theta_3dmat');
if data.use_seasons
    data = rmfield(data, 'dtm_arr');
    data = rmfield(data, 'h_arr');
end

data.theta_summary.illuminated_time_matrix = illuminated_time_matrix;
data.theta_summary.theta_max_matrix = theta_max_matrix;
data.theta_summary.theta_min_matrix = theta_min_matrix;
data.theta_summary.theta_mean_matrix = theta_mean_matrix;
data.theta_summary.source_path = source_path;
data.theta_summary.description = 'Summary of solar angle values (full time series data not included in this file to save space)';
data.theta_summary.summary_created = datetime;
save(target_path, 'data')
fprintf('DONE\n')
end