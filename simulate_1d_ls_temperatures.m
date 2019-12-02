function simulate_1d_ls_temperatures(ls_name)
% SIMULATE_1D_LS_TEMPERATURES simulates landing site temperatures in 1D
% model and saves results in a data file
%
% SIMULATE_1D_LS_TEMPERATURES(ls_name) specifies landing site to simulate

%% Set up 
t_wait = 'seasons';
custom_parameters = struct;

fprintf('Simulating 1D model temperatures for %s\n', ls_name)

%% Load data
source_path = create_static_path(sprintf('crater_environments/%s.mat', ls_name));
crater_data = load(source_path);
crater_data = crater_data.data;
lat_arr = crater_data.lat_arr;
long_arr = crater_data.long_arr;
elevation_matrix = crater_data.elevation_matrix;

if ~isfield(custom_parameters, 'A0') && isfield(crater_data, 'A0')
    custom_parameters.A0 = crater_data.A0;
end

[aspect_matrix, slope_matrix] = convert_height_to_slope(elevation_matrix, lat_arr, long_arr);

Tmax_matrix = NaN(size(elevation_matrix));
Tmin_matrix = NaN(size(elevation_matrix));

[dtm_arr, sub_sun_long_arr, decl_arr] = read_horizons_data();

%% Run simulation
progress_idx = 0;
calculation_t_start_dtm = datetime;
progress_message = '';
for lat_idx = 1:numel(lat_arr)
    lat = lat_arr(lat_idx);
    parfor long_idx = 1:numel(long_arr)
        long = long_arr(long_idx);
        aspect = aspect_matrix(lat_idx, long_idx);
        slope = slope_matrix(lat_idx, long_idx);
        [Tmax, theta_arr] = hayne_1d_model(lat, t_wait, 'max', aspect, slope, custom_parameters, [], long, dtm_arr, sub_sun_long_arr, decl_arr);
        Tmin = hayne_1d_model(lat, t_wait, 'min', aspect, slope, custom_parameters, theta_arr, long, dtm_arr, sub_sun_long_arr, decl_arr);
        Tmax_matrix(lat_idx, long_idx) = Tmax;
        Tmin_matrix(lat_idx, long_idx) = Tmin;
    end
    progress_idx = progress_idx + 1;
    fraction_complete = progress_idx/(numel(lat_arr));%*numel(long_arr));
    finish_dtm = datetime + (datetime - calculation_t_start_dtm)*(1-fraction_complete)/fraction_complete;
    finish_time = datestr(finish_dtm, 'HH:MM:SS dd-mmm');
    simulation_duration = datetime - calculation_t_start_dtm;
    fprintf(repmat('\b', 1, numel(progress_message)))
    progress_message = sprintf('\tCalculation progress: %.1f%%\n\tElapsed time:\t%s\n\tRemaining time:\t%s\n\tPredicted finish: %s\n',...
        100*fraction_complete, string(simulation_duration), string(finish_dtm - datetime), finish_time);
    fprintf('%s', progress_message)
end


data = struct;
data.Tmax_matrix = Tmax_matrix;
data.Tmin_matrix = Tmin_matrix;
data.crater_data = crater_data;
data.t_wait = t_wait;
data.custom_parameters = custom_parameters;
data.description = 'Temperature data for landing site with 1D model';
data.created = datetime;

target_path = create_static_path(sprintf('outputs/ls_1d_temperatures/%s_1d_temperatures.mat', ls_name));
save(target_path, 'data')
end