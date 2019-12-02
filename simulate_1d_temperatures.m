function data = simulate_1d_temperatures(t_wait, T_mode)
% SIMULATE_1D_TEMPERATURES simulates 1D model temperatures for entire lunar
% surface and saves simulated temperatures in a data file.
%
% data = SIMULATE_1D_TEMPERATURES(t_wait, T_mode) simulates 1d model
% temperatues for whole moon with given model wait time and mode, saves and
% returns output data file

if nargin == 0
    simulate_1d_temperatures('match', 'max');
    simulate_1d_temperatures('match', 'min');
    simulate_1d_temperatures(60, 'max');
    simulate_1d_temperatures(60, 'min');
    simulate_1d_temperatures('seasons', 'max');
    simulate_1d_temperatures('seasons', 'min');
    return
end

%% Set up
fprintf('Running 1D model for t_wait=%s, T_mode=%s\n', string(t_wait), string(T_mode))
[height_matrix, lat_arr, long_arr] = read_lola_height_data;
[aspect_matrix, slope_matrix] = convert_height_to_slope(height_matrix, lat_arr, long_arr);

parameters = define_parameters;
albedo_matrix = read_lola_A0_data();
albedo_matrix = scale_data(albedo_matrix, parameters.A0);

if strcmp(t_wait, 'seasons') || strcmp(t_wait, 'match')
    [dtm_arr, sub_sun_long_arr, decl_arr] = read_horizons_data();
    [~, ~, ~, ltim_matrix, dtm_matrix] = read_diviner_temperatures(T_mode);
else
    dtm_arr = [];
    sub_sun_long_arr = [];
    decl_arr = [];
end

Tsimulated_matrix = NaN(numel(lat_arr), numel(long_arr));
simulation_t_start_dtm = datetime;
progress_message = sprintf('\tSimulation progress: 0%%\n\tElapsed time:\t0\n\tRemaining time:\t...\n\tPredicted finish: ...\n');
fprintf('%s', progress_message)
for lat_idx = 1:numel(lat_arr)
    %% Run simulation
    lat = lat_arr(lat_idx);
    parfor long_idx = 1:numel(long_arr)
        long = long_arr(long_idx);
        aspect = aspect_matrix(lat_idx, long_idx);
        slope = slope_matrix(lat_idx, long_idx);
        custom_parameters = struct('A0', albedo_matrix(lat_idx, long_idx));
        if strcmp(t_wait, 'match')
            t_wait_local = dtm_matrix(lat_idx, long_idx);
            T_mode_local = mod(0.5 + ltim_matrix(lat_idx, long_idx)/24,1);
        else
            t_wait_local = t_wait;
            T_mode_local = T_mode;
        end
        Tsimulated_matrix(lat_idx, long_idx) = hayne_1d_model(lat, t_wait_local, T_mode_local, aspect, slope, custom_parameters, [], long, dtm_arr, sub_sun_long_arr, decl_arr);
    end
    
    %% Print progress
    fraction_complete = (lat_idx)/numel(lat_arr);
    finish_dtm = datetime + (datetime - simulation_t_start_dtm)*(1-fraction_complete)/fraction_complete;
    finish_time = datestr(finish_dtm, 'HH:MM:SS dd-mmm');
    simulation_duration = datetime - simulation_t_start_dtm;
    fprintf(repmat('\b', 1, numel(progress_message)))
    progress_message = sprintf('\tSimulation progress: %.1f%%\n\tElapsed time:\t%s\n\tRemaining time:\t%s\n\tPredicted finish: %s\n',...
        100*fraction_complete, string(simulation_duration), string(finish_dtm - datetime), finish_time);
    fprintf('%s', progress_message)
end

%% Save data
data = struct;
data.t_wait = t_wait;
data.T_mode = T_mode;
data.lat_arr = lat_arr;
data.long_arr = long_arr;
data.Tsimulated_matrix = Tsimulated_matrix;
data.description = '1D model temperatures for the entire lunar surface using slopes and scaled LOLA A0';
data.created = datetime;

data_filename = sprintf('outputs/1d_simulated_T/simulated_T_%s_%s.mat', string(t_wait), string(T_mode));
save(data_filename, 'data')

%% Finish
fprintf('\tDone in %s\n', simulation_duration)
end