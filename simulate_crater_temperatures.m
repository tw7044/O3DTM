function data = simulate_crater_temperatures(crater_name, T_mode, t_wait, n_powers, start_power, end_power)
% SIMULATE_CRATER_TEMPERATURES simulates crater temperatures using 3D model
% for different scattering powers and sves results in a data file
%
% data = SIMULATE_CRATER_TEMPERATURES(...) saves crater temperatures and
% returns copy of output data file
%
% SIMULATE_CRATER_TEMPERATURES(crater_name) simulates temperatures for
% specified crater
%
% SIMULATE_CRATER_TEMPERATURES(crater_name, T_mode, t_wait) simulates
% temperatures for specified model arguments
% 
% SIMULATE_CRATER_TEMPERATURES(crater_name, T_mode, t_wait, n_powers)
% specifies number of scattering powers to simulate
% 
% SIMULATE_CRATER_TEMPERATURES(crater_name, T_mode, t_wait, n_powers,
% start_power, end_power) specifies start and end scattering power to
% simulate

%% Set up
if nargin == 0
    crater_name = 'all';
end
if nargin < 2
    T_mode = 'max';
end
if nargin < 3
    t_wait = 'match';
end
if nargin < 4
    n_powers = 11;
end
if nargin < 5
    start_power = 0;
end
if nargin < 6
    end_power = 1;
end
if strcmp(crater_name, 'simulated') || strcmp(crater_name, 'simulated_r')
    crater_name_arr = {'simulated_0n', 'simulated_30n', 'simulated_60n', 'simulated_85n'};
    if strcmp(t_wait, 'match')
        t_wait = 'seasons';
    end
    if strcmp(crater_name, 'all_r')
        % can run simulate_crater_temperatures on 2 computers with 'all'
        % and 'all_r' to run through list from both ends to calculate all
        % craters faster
        crater_name_arr = flip(crater_name_arr);
    end
    for crater_idx = 1:numel(crater_name_arr)
        crater_name = crater_name_arr{crater_idx};
        simulate_crater_temperatures(crater_name, 'max', t_wait, n_powers, start_power, end_power);
        simulate_crater_temperatures(crater_name, 'min', t_wait, n_powers, start_power, end_power);
    end
end
if strcmp(crater_name, 'all') || strcmp(crater_name, 'all_r')
    crater_name_arr = {'blagg', '61n', '88s_compressed', 'bruce', 'erlanger_compressed'};
    if strcmp(crater_name, 'all_r')
        % can run simulate_crater_temperatures on 2 computers with 'all'
        % and 'all_r' to run through list from both ends to calculate all
        % craters faster
        crater_name_arr = flip(crater_name_arr);
    end
    for crater_idx = 1:numel(crater_name_arr)
        crater_name = crater_name_arr{crater_idx};
        simulate_crater_temperatures(crater_name, 'max', t_wait, n_powers, start_power, end_power);
        simulate_crater_temperatures(crater_name, 'min', t_wait, n_powers, start_power, end_power);
    end
    return
end

fprintf('Simulating temperatures: "%s", T_mode=%s, t_wait=%s\n', crater_name, string(T_mode), string(t_wait));
simulation_start_dtm = datetime;
source_path = sprintf('inputs/crater_environments/%s.mat', crater_name);
crater_data = load(source_path);
crater_data = crater_data.data;
scattering_power_arr = linspace(start_power, end_power, n_powers);
Tsimulated_3dmat = NaN(n_powers, numel(crater_data.lat_arr), numel(crater_data.long_arr));


fprintf(['\t', repmat('.',1,n_powers), '\n\t\n']);
parfor power_idx = 1:n_powers
    scattering_power = scattering_power_arr(power_idx);
    T_3dmat = temperature_model_3d(crater_name, T_mode, t_wait, scattering_power, [], struct(), false);
    Tsimulated_3dmat(power_idx, :, :) = T_3dmat(:,:,1);
    fprintf('\b#\n')
end
fprintf('\tSaving data\n')
data = struct;
data.crater_name = crater_data.crater_name;
data.lat_arr = crater_data.lat_arr;
data.long_arr = crater_data.long_arr;
data.ew_matrix = crater_data.ew_matrix;
data.ns_matrix = crater_data.ns_matrix;
data.center_lat = crater_data.center_lat;
data.center_long = crater_data.center_long;
data.elevation_matrix = crater_data.elevation_matrix;
if ~isfield(crater_data, 'simulated_crater') || ~crater_data.simulated_crater
    data.Tmax_diviner_matrix = crater_data.Tmax_matrix;
    data.Tmax_ltim_diviner_matrix = crater_data.Tmax_ltim_matrix;
    data.Tmax_dtm_diviner_matrix = crater_data.Tmax_dtm_matrix;
    data.Tmin_diviner_matrix = crater_data.Tmin_matrix;
    data.Tmin_ltim_diviner_matrix = crater_data.Tmin_ltim_matrix;
    data.Tmin_dtm_diviner_matrix = crater_data.Tmin_dtm_matrix;
    data.ppd = crater_data.ppd;
    data.A0 = crater_data.A0;
else
    data.simulated_crater = true;
end
data.crater_data_created = crater_data.created;
data.scattering_power_arr = scattering_power_arr;
data.Tsimulated_3dmat = Tsimulated_3dmat;
data.T_mode = T_mode;
data.t_wait = t_wait;
data.description = 'Crater temperatures for different scattering power values';
data.created = datetime;

target_name = sprintf('outputs/3d/simulated_T/%s_%s_%s_%ipts.mat', crater_name, string(T_mode), string(t_wait), n_powers);
save(target_name, 'data')

fprintf('\tDone in %s\n', datetime - simulation_start_dtm)
end