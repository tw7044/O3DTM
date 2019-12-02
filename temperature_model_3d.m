function [output, z_arr, crater_data] = temperature_model_3d(crater_name, T_mode, t_wait, scattering_power, num_bounces, custom_parameters, print_progress)
% TEMPERATURE_MODEL_3D 3D thermal model accounting for scattered radiation
% and shadowing effects
%  Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  crater_name = name of crater environment to run simulation for
%  T_mode = controls outpedit analyse_lsut of model:
%  > 'all' => output = struct of all useful results (default)
%  > 'max' => output = maximum lunar surface temperatures
%  > 'min' => output = minimum lunar surface temperatures
%  > number => output = lunar surface temperatures at fraction of lunar day
%  after t_wait
%  t_wait = controls how long model runs for:
%  > number => number of days to run simulation for (to stabilise T) with
%    no seasonal variation
%  > 'seasons' => run with declination data to simulate lunar seasons
%  > 'match' => run to specific time when max/min diviner temperature
%  happened
%  > 'match_ltim' => run with seasonal variation and get max/min
%  temperatures at specific local time that max/min diviner temperature
%  happened at
%  > datetime => run with seasons
%  Optional Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  scattering_power = emissivity scattering power (default = 0)
%  num_bounces = max. number of bounces for scattered radiation (default =
%  inf)
%  custom_parameters = struct of custom parameters for code, A0, a, b, H
%  and emissivity can be matrices (e.g. albedo varying across map)
%  print_progress = time step interval to print progress, 0 supresses all
%  output, 1 prints progress for every time step, a larger integer results
%  in (slightly) faster execution (default = 360)
%
%  Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Return max temperatures ocurring in 'blagg' crater accounting for
%  seasonal variation:
%    temperature_model_3d('blagg', 'max', 'seasons')
%  Return temperatures occurring in 'bruce' crater at the datetime that the
%  minimum Diviner measurements were made:
%    temperature_model_3d('bruce', 'min', 'match')
%  Return temperature ocurring in the simulated crater at the equator after
%  10 lunar days of simulation:
%    temperature_model_3d('simulated_0n', 0, 10)
%
% By Oliver King, 2018

narginchk(1,inf)
simulation_start_dtm = datetime;

%% Process inputs
if nargin < 2 || isempty(T_mode)
    T_mode = 'all';
end
if nargin < 3 || isempty(t_wait)
    t_wait = 60;
end
if nargin < 4 || isempty(scattering_power)
    scattering_power = 0; % default to isotropic emission
end
if nargin < 5 || isempty(num_bounces)
    num_bounces = inf;
end
if nargin < 6 || isempty(custom_parameters)
    custom_parameters = struct();
end
if nargin < 7 || isempty(print_progress)
    print_progress = 360;
end

if ischar(T_mode) && strcmp(T_mode, 'all')
    t_wait = 'seasons';
end

if print_progress
    T_mode_str = string(T_mode);
    t_wait_str = string(t_wait);
    fprintf('3D model: "%s", T_mode=%s, t_wait=%s, scattering_power=%g, num_bounces=%g\n', crater_name, T_mode_str, t_wait_str, scattering_power, num_bounces)
end

%% Load crater data
source_path = create_static_path(sprintf('crater_environments/%s.mat', crater_name));
try
    crater_data = load(source_path);
catch ME
    if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
        crater_data = generate_crater_environment(crater_name);
    else
        error(ME.message)
    end
end
crater_data = crater_data.data;
lat_arr = crater_data.lat_arr;
long_arr = crater_data.long_arr;
ew_matrix = crater_data.ew_matrix;
ns_matrix = crater_data.ns_matrix;
center_lat = crater_data.center_lat;
center_long = crater_data.center_long;
elevation_matrix = crater_data.elevation_matrix;
if ~isfield(custom_parameters, 'A0') && isfield(crater_data, 'A0')
    custom_parameters.A0 = crater_data.A0;
end


%% Process info from inputs
% Get information about terrain and linear distances etc.
ref_sphere = referenceSphere('moon');
r_moon = ref_sphere.Radius;
elevation_matrix = elevation_matrix - r_moon;

[aspect_matrix, slope_matrix] = convert_height_to_slope(elevation_matrix, lat_arr, long_arr);

% load appropriate angle information
if (ischar(t_wait) && (strcmp(t_wait, 'seasons') || strcmp(t_wait, 'match'))) || strcmp(t_wait, 'match_ltim') || isdatetime(t_wait)
    use_seasons = true;
    file_season_str = 'seasons';
else
    use_seasons = false;
    file_season_str = 'no_seasons';
end
source_path = create_static_path(sprintf('crater_environments/ray_tracing/%s/%s_ray_tracing.mat', file_season_str, crater_name));
% re-use theta_3dmat variable name to save memory (as we don't need to keep
% the struct after extracting the useful information from it)
theta_3dmat = load(source_path);
theta_3dmat = theta_3dmat.data;
vis_4dmat = theta_3dmat.vis_4dmat;
ray_tracing_datetime = theta_3dmat.created;
if use_seasons
    time_series_dtm_arr = theta_3dmat.dtm_arr;
    h_arr = theta_3dmat.h_arr;
    custom_parameters.dt = theta_3dmat.dt;
end
theta_3dmat = theta_3dmat.theta_3dmat;

%% Load parameters
parameters = define_parameters(custom_parameters);
% load parameters from define_parameters and convert to individual
% variables to increase code legibility and performance (e.g. using
% parameters.P is slower than using P)
P = parameters.P;
rho_s = parameters.rho_s;
rho_d = parameters.rho_d;
H_matrix = parameters.H;
stefans_constant = parameters.stefans_constant;
emissivity_matrix = parameters.emissivity;
S = parameters.S;
R_AU = parameters.R_AU;
Q = parameters.Q;
A0_matrix = parameters.A0;
a_matrix = parameters.a;
b_matrix = parameters.b;
Chi = parameters.Chi;
Ks = parameters.Ks;
Kd = parameters.Kd;
c0 = parameters.c0;
c1 = parameters.c1;
c2 = parameters.c2;
c3 = parameters.c3;
c4 = parameters.c4;
m = parameters.m;
n = parameters.n;
num_skin_depths = parameters.num_skin_depths;
max_depth = parameters.max_depth;
live_graph = parameters.live_graph;
live_graph_plot_interval = parameters.live_graph_plot_interval;
T_min = parameters.T_min;
T_max = parameters.T_max;
T_bottom_limit = parameters.T_bottom_limit;
surface_bc_test_difference = parameters.surface_bc_test_difference;
surface_bc_break_counter = parameters.surface_bc_break_counter;
dt = parameters.dt;
diviner_start_dtm = parameters.diviner_start_dtm;
diviner_end_dtm = parameters.diviner_end_dtm;
time_series_interval = parameters.time_series_interval;
initial_depth_t_wait = parameters.initial_depth_t_wait;
grow_depth_t_wait = parameters.grow_depth_t_wait;
initial_num_skin_depths = parameters.initial_num_skin_depths;


% generate matrices if needed
if numel(A0_matrix) == 1
    A0_matrix = A0_matrix*ones(size(elevation_matrix));
end
if numel(a_matrix) == 1
    a_matrix = a_matrix*ones(size(elevation_matrix));
end
if numel(b_matrix) == 1
    b_matrix = b_matrix*ones(size(elevation_matrix));
end
if numel(emissivity_matrix) == 1
    emissivity_matrix = emissivity_matrix*ones(size(elevation_matrix));
end
if numel(H_matrix) == 1
    H_matrix = H_matrix*ones(size(elevation_matrix));
end

% pre-calculate constant parameters
S_OVER_R_AU2 = S/R_AU^2;
thermal_emission_matrix = emissivity_matrix*stefans_constant;
Chi_OVER_350_POWER_3 = Chi*(1/350)^3;
a_OVER_pi_OVER_4_POWER_3_matrix = a_matrix*(1/(pi/4))^3;
b_OVER_pi_OVER_2_POWER_8_matrix = b_matrix*(1/(pi/2))^8;


%% Generate grid
% Use smallest expected skin depth, zs_min to define thickness of layers
% and use largest expected skin depth, zs_max to define total number of
% layers to use

    function [K_value] = K_function(T, rho, par)
        Kc_value = par.Kd - (par.Kd - par.Ks)*(par.rho_d - rho)/(par.rho_d - par.rho_s);
        K_value = Kc_value*(1 + par.Chi*(T/350)^3);
    end
    function [c_p_value] = c_p_function(T, par)
        c_p_value = par.c0 + par.c1*T + par.c2*T^2 + par.c3*T^3 + par.c4*T^4;
    end

zs_max = sqrt((K_function(T_max, rho_d, parameters)/(rho_s*c_p_function(T_max, parameters)))*P/pi);
zs_min = sqrt((K_function(T_min, rho_s, parameters)/(rho_s*c_p_function(T_min, parameters)))*P/pi);
layer_idx = 0;
dz = zs_min/m;
z = 0;
max_depth_bool = true;
while (z < zs_max*(num_skin_depths + 1)) || max_depth_bool
    layer_idx = layer_idx + 1;
    z_arr(layer_idx) = z;
    dz_arr(layer_idx) = dz;
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            if H_matrix(lat_idx, long_idx) == 0 && z == 0
                rho_matrix(lat_idx, long_idx, layer_idx) = rho_s;
            else
                rho_matrix(lat_idx, long_idx, layer_idx) = rho_d - (rho_d - rho_s)*exp(-z/H_matrix(lat_idx, long_idx));
            end
        end
    end
    dz = dz*(1+1/n);
    if z > max_depth
        max_depth_bool = false; % ensure have both at least max depth & sufficient skin depths
    end
    z = z + dz;
end

num_layers = layer_idx;

initial_depth = zs_max*initial_num_skin_depths;
bottom_layer_idx = round(interp1(z_arr, 1:numel(z_arr), initial_depth));
depth_update_wait_t = initial_depth_t_wait;
depth_update_t_interval = grow_depth_t_wait/(num_layers - bottom_layer_idx - 1);


%% Calculate simulation duration
if strcmp(t_wait, 'match')
    % prepare recording of temperature for each pixel at time specified by
    % measurement ltim and dtm where simulation runs to first ltim after
    % the specified dtm
    switch T_mode
        case 'max'
            px_ltim_matrix = crater_data.Tmax_ltim_matrix;
            px_dtm_matrix = crater_data.Tmax_dtm_matrix;
        case 'min'
            px_ltim_matrix = crater_data.Tmin_ltim_matrix;
            px_dtm_matrix = crater_data.Tmin_dtm_matrix;
    end
    
    % convert ltim to hour angle
    px_h_matrix = (px_ltim_matrix - 12)*(360/24);
    px_h_matrix = mod(px_h_matrix + 180, 360) - 180; % ensure consistent range
    
    start_dtm = time_series_dtm_arr(1);
    update_t_idx_matrix = NaN(size(elevation_matrix));
    
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            % adjust hour angle array for actual local time depending on
            % long of current pixel
            long_offset = long_arr(long_idx) - center_long;
            local_h_arr = h_arr + long_offset;
            local_h_arr = mod(local_h_arr + 180, 360) - 180;
            
            measure_start_dtm = px_dtm_matrix(lat_idx, long_idx);
            end_dtm = px_dtm_matrix(lat_idx, long_idx) + seconds(P*1.5);
            local_dtm_arr = start_dtm:seconds(dt):end_dtm;
            
            h_end = px_h_matrix(lat_idx, long_idx);

            idx = round(interp1(local_dtm_arr, 1:numel(local_dtm_arr), measure_start_dtm));
            h_end_big = h_end > max(h_arr);
            h_end_small = h_end < min(h_arr);
            for idx = idx:numel(h_arr)
                if (h_arr(idx-1) < h_end) && (h_arr(idx) >= h_end)
                    % on rising slope
                    break
                elseif (h_end_big) && (h_arr(idx-1) > h_arr(idx))
                    % h_end = 180
                    break
                elseif (h_end_small) && (h_arr(idx-1) > h_arr(idx))
                    % h_end = -180
                    break
                end
            end
            update_t_idx_matrix(lat_idx, long_idx) = idx;
        end
    end    
    end_dtm = start_dtm + seconds((max(update_t_idx_matrix(:)))*dt);
    t_limit = seconds(end_dtm - start_dtm);
elseif use_seasons
    % fixed end time
    start_dtm = time_series_dtm_arr(1);
    
    if numel(T_mode) > 1
        t_wait_matrix = t_wait;
    end
    
    if isdatetime(t_wait) && isnumeric(T_mode)
        % run to date then run to ltim
        measure_start_dtm = t_wait;
        end_dtm = t_wait + seconds(P*1.5);
        t_arr = start_dtm:seconds(dt):end_dtm;
        
        h_end = T_mode*360/24;
        h_end = mod(h_end + 180, 360) - 180;
        
        idx = round(interp1(t_arr, 1:numel(t_arr), measure_start_dtm));
        h_end_big = h_end > max(h_arr);
        h_end_small = h_end < min(h_arr);
        for idx = idx:numel(h_arr)
            if (h_arr(idx-1) < h_end) && (h_arr(idx) >= h_end)
                % on rising slope
                break
            elseif (h_end_big) && (h_arr(idx-1) > h_arr(idx))
                % h_end = 180
                break
            elseif (h_end_small) && (h_arr(idx-1) > h_arr(idx))
                % h_end = -180
                break
            end
        end
        end_dtm = t_arr(idx);
        t_wait = seconds(t_wait - start_dtm)/P;
    elseif isdatetime(T_mode)
        % run to date
        end_dtm = T_mode;
        t_wait = 0;
    else
        % run over whole period (max/min)
        end_dtm = diviner_end_dtm;
        t_wait = seconds(diviner_start_dtm - start_dtm)/P;
    end
    t_limit = seconds(end_dtm - start_dtm);
    
    
    if strcmp(t_wait, 'match_ltim')
        % get max/min temperature from seasonal cycle at specified ltim
        switch T_mode
            case 'max'
                px_ltim_matrix = crater_data.Tmax_ltim_matrix;
                px_dtm_matrix = crater_data.Tmax_dtm_matrix;
            case 'min'
                px_ltim_matrix = crater_data.Tmin_ltim_matrix;
                px_dtm_matrix = crater_data.Tmin_dtm_matrix;
        end
        
        % convert ltim to hour angle
        px_measure_h_matrix = (px_ltim_matrix - 12)*(360/24);
        px_measure_h_matrix = mod(px_measure_h_matrix + 180, 360) - 180; % ensure consistent range
        
        long_step_h_matrix = NaN(numel(long_arr), numel(h_arr));
        
        for long_idx = 1:numel(long_arr)
            % adjust hour angle array for actual local time depending on
            % long of current pixel
            long_offset = long_arr(long_idx) - center_long;
            local_h_arr = h_arr + long_offset;
            
            local_h_arr = mod(local_h_arr + 180, 360) - 180;
            
            long_step_h_matrix(long_idx, :) = local_h_arr;
        end
    end
else
    % non-seasonal simulation
    if ischar(T_mode)
        t_limit = P*(t_wait + 1) - dt;
    else
        t_limit = P*(t_wait + T_mode);
    end
end


%% Pre-calculate information about each pixel
if print_progress
    fprintf('\tCalculating flux coefficients')
    fprintf([repmat('.',1,numel(lat_arr)-numel('Calculating flux coefficients')), '\n\t\n']);
end

absorbed_vis_coeff_4dmat = zeros(numel(lat_arr), numel(long_arr), numel(lat_arr), numel(long_arr));
scattered_vis_coeff_4dmat = zeros(size(absorbed_vis_coeff_4dmat));
emitted_absorbed_ir_coeff_4dmat = zeros(size(absorbed_vis_coeff_4dmat));
emitted_scattered_ir_coeff_4dmat = zeros(size(absorbed_vis_coeff_4dmat));
scattered_absorbed_ir_coeff_4dmat = zeros(size(absorbed_vis_coeff_4dmat));
scattered_scattered_ir_coeff_4dmat = zeros(size(absorbed_vis_coeff_4dmat));

% calculate normalisations to ensure energy conservation for non-isotropic
% scattering functions
syms theta
ir_normalisation = 1/double(int(sin(theta)*(cos(theta))^scattering_power, theta, 0, pi/2)); % = scattering_power + 1
vis_normalisation_matrix = 1./double(int(sin(theta).*(A0_matrix + a_matrix.*(theta/(pi/4))^3 + b_matrix.*(theta/(pi/2))^8), theta, 0, pi/2)); % = ~3 for default parameters
clear theta

T_3dmat = NaN(numel(lat_arr), numel(long_arr), num_layers);
T_new_3dmat = NaN(size(T_3dmat));

for lat1_idx = 1:numel(lat_arr)
    lat1 = lat_arr(lat1_idx);

    %% Prepare for parfor iteration
    % Set up local copies so parfor works (if using parfor version of code)
    T_3dmat_local = T_3dmat(lat1_idx, :, :);
    absorbed_vis_coeff_4dmat_local = absorbed_vis_coeff_4dmat(lat1_idx, :, :, :);
    scattered_vis_coeff_4dmat_local = scattered_vis_coeff_4dmat(lat1_idx, :, :, :);
    
    emitted_absorbed_ir_coeff_4dmat_local = emitted_absorbed_ir_coeff_4dmat(lat1_idx, :, :, :);
    emitted_scattered_ir_coeff_4dmat_local = emitted_scattered_ir_coeff_4dmat(lat1_idx, :, :, :);
    scattered_absorbed_ir_coeff_4dmat_local = scattered_absorbed_ir_coeff_4dmat(lat1_idx, :, :, :);
    scattered_scattered_ir_coeff_4dmat_local = scattered_scattered_ir_coeff_4dmat(lat1_idx, :, :, :);
    
    for long1_idx = 1:numel(long_arr)
        long1 = long_arr(long1_idx);
        
        %% Get information about current location
        slope1 = slope_matrix(lat1_idx, long1_idx);
        aspect1 = aspect_matrix(lat1_idx, long1_idx);
        h1 = elevation_matrix(lat1_idx, long1_idx);
                
        %% Initialise temperature profile
        % Assume surface in eqilibrium with F_solar and black body
        % emission, and deep temperature drops off to T0/sqrt(2). Assumes
        % the temperature at the start is due to the most direct solar
        % radiation on the surface during the simulation period.
        theta = double(nanmin(theta_3dmat(lat1_idx, long1_idx, :)));
        albedo = A0_matrix(lat1_idx, long1_idx) + a_OVER_pi_OVER_4_POWER_3_matrix(lat1_idx, long1_idx)*theta*theta*theta + b_OVER_pi_OVER_2_POWER_8_matrix(lat1_idx, long1_idx)*theta*theta*theta*theta*theta*theta*theta*theta;
        T0 = (((1-albedo)/(emissivity_matrix(lat1_idx, long1_idx)*stefans_constant))*S/R_AU^2*cos(theta))^0.25;
        if theta > pi/2 || T0 < T_bottom_limit || ~isreal(T0) || isnan(T0)
            T0 = T_bottom_limit;
        end
        TN = T0/sqrt(2);
        if TN < T_bottom_limit
            % set roughly constant initial temperature profile for very
            % cold regions
            TN = T_bottom_limit;
        end
        for layer_idx = 1:num_layers
            z = z_arr(layer_idx);
            if H_matrix(lat1_idx, long1_idx) == 0 && z == 0
                T_3dmat_local(1, long1_idx, layer_idx) = T0;
            else
                T_3dmat_local(1, long1_idx, layer_idx) = TN - (TN - T0)*exp(-z/H_matrix(lat1_idx, long1_idx));
            end
        end
        
        %% Calculate flux coefficients
        % Calculate which pixels are visible from the current pixel and
        % calculate values to be used in flux coeff calculation
        vis_matrix = squeeze(vis_4dmat(lat1_idx, long1_idx, :, :));
        ew_dist = ew_matrix(lat1_idx, 2) - ew_matrix(lat1_idx, 1);
        ns_dist = ns_matrix(1, long1_idx) - ns_matrix(2, long1_idx);
        emission_area = abs(ew_dist*ns_dist/cosd(slope_matrix(lat1_idx, long1_idx)));
        
        % Calculate flux coefficients for each pixel visible from current
        % pixel (i.e. radiation emitted from px1 which is absorbed by px2
        % and which is scattered from px2)
        for lat2_idx = 1:numel(lat_arr)
            lat2 = lat_arr(lat2_idx);
            for long2_idx = 1:numel(long_arr)
                long2 = long_arr(long2_idx);
                if lat1 == lat2 && long1 == long2
                    % Values not needed for current pixel
                    continue
                elseif vis_matrix(lat2_idx, long2_idx) == 0
                    % Values not needed for pixels not visible from current
                    % pixel
                    continue
                else
                    h2 = elevation_matrix(lat2_idx, long2_idx);
                    [az12,elev12,dist] = geodetic2aer(lat1,long1,h1,lat2,long2,h2,ref_sphere);
                    [az21,elev21,~] = geodetic2aer(lat2,long2,h2,lat1,long1,h1,ref_sphere);
                    % calculate flux coefficients          
                    slope2 = slope_matrix(lat2_idx, long2_idx);
                    aspect2 = aspect_matrix(lat2_idx, long2_idx);
                    
                    elev12 = 90 - abs(elev12);
                    elev21 = 90 - abs(elev21);

                    cos_theta1 = cosd(elev12)*cosd(slope1) + sind(elev12)*sind(slope1)*cosd(az12-aspect1);
                    cos_theta2 = cosd(elev21)*cosd(slope2) + sind(elev21)*sind(slope2)*cosd(az21-aspect2);
                    
                    theta1 = acos(cos_theta1);
                    theta2 = acos(cos_theta2);
                    
                    %% SCATTERING FUNCTION #############################################################################################
                    % geometric_flux_coeff is fraction of radiation emitted
                    % from px1 which can be absorbed by px2 (assuming
                    % perfect zero albedo)
                    % absorbed_flux_coeff is fraction of radiation emitted
                    % from px1 which is absorbed by px2 (accounting for
                    % albedo)
                    % scattered_flux_coeff is fraction of radiation from
                    % px1 which is incident on px2 but is re-scattered due
                    % to albedo
                    
                    %% Geometric
                    geometric_flux_coeff = cos_theta1*cos_theta2*emission_area/(4*pi*dist^2);
                    
                    if geometric_flux_coeff < 0 || geometric_flux_coeff > 1 || cos_theta1 < 0 || cos_theta2 < 0 || ~isreal(geometric_flux_coeff)
                        % run check to ensure physical coefficient
                        geometric_flux_coeff = 0;
                    end

                    %% Visible
                    albedo1 = A0_matrix(lat1_idx, long1_idx) + a_OVER_pi_OVER_4_POWER_3_matrix(lat1_idx, long1_idx)*theta1*theta1*theta1 + b_OVER_pi_OVER_2_POWER_8_matrix(lat1_idx, long1_idx)*theta1*theta1*theta1*theta1*theta1*theta1*theta1*theta1;
                    albedo2 = A0_matrix(lat2_idx, long2_idx) + a_OVER_pi_OVER_4_POWER_3_matrix(lat2_idx, long2_idx)*theta2*theta2*theta2 + b_OVER_pi_OVER_2_POWER_8_matrix(lat2_idx, long2_idx)*theta2*theta2*theta2*theta2*theta2*theta2*theta2*theta2;
                    
                    absorbed_vis_coeff = vis_normalisation_matrix(lat1_idx, long1_idx)*albedo1*geometric_flux_coeff*(1-albedo2); % scattered from 1, absorbed by 2
                    scattered_vis_coeff = vis_normalisation_matrix(lat1_idx, long1_idx)*albedo1*geometric_flux_coeff*(albedo2); % scattered from 1, scattered by 2
                    
                    absorbed_vis_coeff_4dmat_local(1, long1_idx, lat2_idx, long2_idx) = absorbed_vis_coeff;                    
                    scattered_vis_coeff_4dmat_local(1, long1_idx, lat2_idx, long2_idx) = scattered_vis_coeff;
                    
                    %% IR
                    emissivity1 = emissivity_matrix(lat1_idx, long1_idx)*cos_theta1^scattering_power;
                    emissivity2 = emissivity_matrix(lat2_idx, long2_idx)*cos_theta2^scattering_power;
                    
                    emitted_absorbed_ir_coeff = ir_normalisation*(emissivity1)*geometric_flux_coeff*(emissivity2); % emitted from 1, absorbed by 2
                    emitted_scattered_ir_coeff = ir_normalisation*(emissivity1)*geometric_flux_coeff*(1-emissivity2); % emitted from 1, scattered by 2
                    scattered_absorbed_ir_coeff = (1-emissivity1)*geometric_flux_coeff*(emissivity2); % scattered from 1, absorbed by 2
                    scattered_scattered_ir_coeff = (1-emissivity1)*geometric_flux_coeff*(1-emissivity2); % scattered from 1, scattered by 2
                    
                    emitted_absorbed_ir_coeff_4dmat_local(1, long1_idx, lat2_idx, long2_idx) = emitted_absorbed_ir_coeff;
                    emitted_scattered_ir_coeff_4dmat_local(1, long1_idx, lat2_idx, long2_idx) = emitted_scattered_ir_coeff;
                    scattered_absorbed_ir_coeff_4dmat_local(1, long1_idx, lat2_idx, long2_idx) = scattered_absorbed_ir_coeff;
                    scattered_scattered_ir_coeff_4dmat_local(1, long1_idx, lat2_idx, long2_idx) = scattered_scattered_ir_coeff;

                    %% #################################################################################################################
                end
            end
        end
    end
    
    %% Finish parfor iteration
    % Finish parfor iteration by saving local copies
    T_3dmat(lat1_idx, :, :) = T_3dmat_local;
    absorbed_vis_coeff_4dmat(lat1_idx, :, :, :) = absorbed_vis_coeff_4dmat_local;
    scattered_vis_coeff_4dmat(lat1_idx, :, :, :) = scattered_vis_coeff_4dmat_local;
    emitted_absorbed_ir_coeff_4dmat(lat1_idx, :, :, :) = emitted_absorbed_ir_coeff_4dmat_local;
    emitted_scattered_ir_coeff_4dmat(lat1_idx, :, :, :) = emitted_scattered_ir_coeff_4dmat_local;
    scattered_absorbed_ir_coeff_4dmat(lat1_idx, :, :, :) = scattered_absorbed_ir_coeff_4dmat_local;
    scattered_scattered_ir_coeff_4dmat(lat1_idx, :, :, :) = scattered_scattered_ir_coeff_4dmat_local;
    % Update progress bar
    if print_progress
        fprintf('\b#\n')
    end
end


%% Perform multiple scattering calculations
% Re-index as 2d matrix
num_px = numel(lat_arr)*numel(long_arr);
absorbed_vis_flux_matrix = reshape(absorbed_vis_coeff_4dmat, num_px, num_px);
scattered_vis_flux_matrix = reshape(scattered_vis_coeff_4dmat, size(absorbed_vis_flux_matrix));
emitted_absorbed_ir_flux_matrix = reshape(emitted_absorbed_ir_coeff_4dmat, size(absorbed_vis_flux_matrix));
emitted_scattered_ir_flux_matrix = reshape(emitted_scattered_ir_coeff_4dmat, size(absorbed_vis_flux_matrix));
scattered_absorbed_ir_flux_matrix = reshape(scattered_absorbed_ir_coeff_4dmat, size(absorbed_vis_flux_matrix));
scattered_scattered_ir_flux_matrix = reshape(scattered_scattered_ir_coeff_4dmat, size(absorbed_vis_flux_matrix));
A0_nx1 = reshape(A0_matrix, [], 1);
a_OVER_pi_OVER_4_POWER_3_nx1 = reshape(a_OVER_pi_OVER_4_POWER_3_matrix, [], 1);
b_OVER_pi_OVER_2_POWER_8_nx1 = reshape(b_OVER_pi_OVER_2_POWER_8_matrix, [], 1);

% save memory for large craters
clear absorbed_vis_coeff_4dmat scattered_vis_coeff_4dmat
clear emitted_absorbed_ir_coeff_4dmat emitted_scattered_ir_coeff_4dmat
clear emitted_absorbed_ir_coeff_4dmat emitted_scattered_ir_coeff_4dmat

if num_bounces < 0
    % Re-absorb radiation into surface
    num_bounces = -num_bounces;
    emitted_absorbed_ir_flux_matrix = emitted_absorbed_ir_flux_matrix/ir_normalisation;
    emitted_scattered_ir_flux_matrix = emitted_scattered_ir_flux_matrix/ir_normalisation;
    emitted_absorbed_ir_flux_matrix = emitted_absorbed_ir_flux_matrix + eye(size(emitted_absorbed_ir_flux_matrix))*(ir_normalisation-1)/ir_normalisation;    
end

if num_bounces == 0
    vis_flux_matrix = zeros(size(absorbed_vis_flux_matrix));
    ir_flux_matrix = zeros(size(emitted_absorbed_ir_flux_matrix));
elseif num_bounces == 1
    vis_flux_matrix = absorbed_vis_flux_matrix;
    ir_flux_matrix = emitted_absorbed_ir_flux_matrix;
else
    if isinf(num_bounces)
        multiple_vis_scatter_matrix = (eye(size(scattered_vis_flux_matrix)) - scattered_vis_flux_matrix)^-1; % geometric series sum to infinity
        multiple_ir_scatter_matrix = (eye(size(scattered_scattered_ir_flux_matrix)) - scattered_scattered_ir_flux_matrix)^-1; % geometric series sum to infinity
    else
        multiple_vis_scatter_matrix = (eye(size(scattered_vis_flux_matrix)) - scattered_vis_flux_matrix^(num_bounces-1))*(eye(size(scattered_vis_flux_matrix)) - scattered_vis_flux_matrix)^-1;
        multiple_ir_scatter_matrix = (eye(size(scattered_scattered_ir_flux_matrix)) - scattered_scattered_ir_flux_matrix^(num_bounces-1))*(eye(size(scattered_scattered_ir_flux_matrix)) - scattered_scattered_ir_flux_matrix)^-1;
    end
    vis_flux_matrix = absorbed_vis_flux_matrix + scattered_vis_flux_matrix*multiple_vis_scatter_matrix*absorbed_vis_flux_matrix;
    ir_flux_matrix = emitted_absorbed_ir_flux_matrix + emitted_scattered_ir_flux_matrix*multiple_ir_scatter_matrix*scattered_absorbed_ir_flux_matrix;
end
pixel_idx_matrix = reshape(1:num_px, size(elevation_matrix)); % allows switching between lat/long indices and pixel index
% deal with floating point errors etc. to ensure flux coefficents never
% negative
vis_flux_matrix(vis_flux_matrix < 0) = 0;
ir_flux_matrix(ir_flux_matrix < 0) = 0;
clear absorbed_vis_flux_matrix scattered_vis_flux_matrix
clear emitted_absorbed_ir_flux_matrix emitted_scattered_ir_flux_matrix
clear emitted_absorbed_ir_flux_matrix emitted_scattered_ir_flux_matrix

%% Prepare for loop
% Prepare for loop by pre-calculating constant values which do not need to
% be calculated every iteration (therefore increasing performance)
Kc_3dmat = NaN(numel(lat_arr), numel(long_arr), num_layers);
for lat_idx = 1:numel(lat_arr)
    for long_idx = 1:numel(long_arr)
        for layer_idx = 1:num_layers
            Kc_3dmat(lat_idx, long_idx, layer_idx) = Kd - (Kd - Ks)*(rho_d - rho_matrix(lat_idx, long_idx, layer_idx))/(rho_d - rho_s);
        end
        B_surface_matrix(lat_idx, long_idx) = Kc_3dmat(lat_idx, long_idx, 1)*Chi_OVER_350_POWER_3;
    end
end

p_arr = NaN(1,num_layers);
q_arr = NaN(1,num_layers);
for layer_idx = 2:num_layers-1
    d3z = dz_arr(layer_idx)*dz_arr(layer_idx-1)*(dz_arr(layer_idx) + dz_arr(layer_idx-1));
    p = 2*dz_arr(layer_idx)/d3z;
    q = 2*dz_arr(layer_idx-1)/d3z;
    p_arr(layer_idx) = p;
    q_arr(layer_idx) = q;
end


%% Perform temperature simulation
if print_progress
    fprintf('\tRunning simulation\n')
end

T_new_arr = NaN(1,num_layers);
t_steps = P/dt;
t_arr = 0:dt:t_limit;
dh = 360/t_steps;

if ischar(T_mode) && strcmp(T_mode, 'all')
    Tmax_3dmat = NaN(size(T_3dmat));
    Tmax_dtm_3dmat = NaN(size(T_3dmat));
    Tmin_3dmat = NaN(size(T_3dmat));
    Tmin_dtm_3dmat = NaN(size(T_3dmat));
    
    time_series_start = t_limit - 60*60*24*365.25; % 1 (Julian) year
    time_series_dtm_arr = start_dtm + seconds(t_arr((t_arr >= time_series_start) & (mod(1:numel(t_arr), time_series_interval) == 0)));
    T_time_series_4dmat = NaN([size(T_3dmat), numel(time_series_dtm_arr)], 'single');
    time_series_idx = 1;
else
    T_output_3dmat = NaN(size(T_3dmat));
end

if live_graph
    live_graph_idx = 1;
end
simulation_t_start_dtm = datetime;
progress_message = '';
for t_idx = 1:numel(t_arr)
    t = t_arr(t_idx);
    if use_seasons
        t_idx_local = t_idx;
    else
        h = 360*mod(t,P)/P;
        t_idx_local = floor(h/dh + 1);
    end
    
    %% Increase simulation depth
    if bottom_layer_idx < num_layers && t > depth_update_wait_t
        % extend simulation by 1 layer, using bottom BC to define new
        % temperature
        depth_update_wait_t = depth_update_wait_t + depth_update_t_interval;
        bottom_layer_idx = bottom_layer_idx + 1;
        T_above_matrix = T_3dmat(:,:,bottom_layer_idx-1);
        T_3dmat(:,:,bottom_layer_idx) = T_above_matrix + (dz_arr(bottom_layer_idx-1)*Q)./(Kc_3dmat(:,:,bottom_layer_idx-1).*(1 + Chi_OVER_350_POWER_3.*T_above_matrix.*T_above_matrix.*T_above_matrix));
    end
    
    %% Process matrix values
    theta_nx1 = reshape(double(theta_3dmat(:,:,t_idx_local)), [], 1); % convert back up from single precision
    solar_albedo_nx1 = A0_nx1 + a_OVER_pi_OVER_4_POWER_3_nx1.*theta_nx1.*theta_nx1.*theta_nx1 + b_OVER_pi_OVER_2_POWER_8_nx1.*theta_nx1.*theta_nx1.*theta_nx1.*theta_nx1.*theta_nx1.*theta_nx1.*theta_nx1.*theta_nx1;
    Qsun_nx1 = (1 - solar_albedo_nx1).*S_OVER_R_AU2.*cos(theta_nx1); % direct solar flux absorbed by each px
    Qsun_scattered_nx1 = solar_albedo_nx1.*S_OVER_R_AU2.*cos(theta_nx1);
    Qsun_nx1(isnan(Qsun_nx1)) = 0; % deal with NaN values represention no solar illumination
    Qsun_scattered_nx1(isnan(Qsun_scattered_nx1)) = 0;
    sigmaT4_nx1 = reshape(stefans_constant*T_3dmat(:,:,1).^4, [], 1);
    Qabsorbed_vis_nx1 = vis_flux_matrix*Qsun_scattered_nx1; % scattered vis flux absorbed by each px
    Qabsorbed_ir_nx1 = ir_flux_matrix*sigmaT4_nx1; % scattered ir flux absorbed by each px
    Q_arr = Qabsorbed_ir_nx1 + Qabsorbed_vis_nx1 + Qsun_nx1; % total flux absorbed by each px
    
    %% Live graph plotting
    % Plot live graph if enabled (generally for debugging purposes as this
    % will decrease execution speed my many orders of magnitude)
    if live_graph
%         if mod(live_graph_idx, live_graph_plot_interval) == 0 || t == 0 || t + dt > t_limit
%             clf
%             hold on
%             lat_idx = 3;
%             long_idx = 3;
%             plot(squeeze(Tmin_3dmat(lat_idx, long_idx, :)), -z_arr)
%             plot(squeeze(Tmax_3dmat(lat_idx, long_idx, :)), -z_arr)
%             plot(squeeze(T_3dmat(lat_idx, long_idx, :)), -z_arr, 'k')
%             xlabel('Temperature (K)')
%             ylabel('Depth (m)')
%             title(sprintf('Layer temperatures, simulation duration = %0.2f lunar days', t/P))
%             grid on
%             drawnow
%         end
        
        if mod(live_graph_idx, live_graph_plot_interval) == 0 || t == 0 || t + dt >= t_limit
            if ~exist('T_output_3dmat', 'var') || numel(T_output_3dmat) - sum(isnan(T_output_3dmat(:))) == 0
                plot_matrix = T_3dmat(:,:,1);
                plot_title = sprintf('Live temperature (%0.3f)', t/P);
            else
                plot_matrix = T_output_3dmat(:,:,1);
                plot_title = sprintf('Output temperature (%0.3f)', t/P);
            end
            plot_matrix = plot_matrix - crater_data.Tmax_matrix;
            % plot_matrix = reshape((Qabsorbed_ir_nx1 + Qabsorbed_vis_nx1)./Q_arr, size(elevation_matrix));
            plot_3d_surface(ew_matrix, ns_matrix, elevation_matrix, plot_matrix, plot_title);
            % divergent_colormap(100);
            step_colormap(spectral, 1, 0, true);
            drawnow
        end
        live_graph_idx = live_graph_idx + 1;
    end
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            T_arr = T_3dmat(lat_idx, long_idx, :);
            pixel_idx = pixel_idx_matrix(lat_idx, long_idx);
            %% Surface BC
            break_iter = 0;
            difference = surface_bc_test_difference + 1;
            T = T_arr(1);
            Q_surface = Q_arr(pixel_idx);
            % Root finding
            while difference > surface_bc_test_difference
                K_value = Kc_3dmat(lat_idx, long_idx, 1)*(1 + Chi_OVER_350_POWER_3*T*T*T);
                func = thermal_emission_matrix(lat_idx, long_idx)*T*T*T*T - Q_surface - K_value*((-3*T + 4*T_arr(2) - T_arr(3))/(2*dz_arr(1)));
                diff = 4*thermal_emission_matrix(lat_idx, long_idx)*T*T*T - 3*B_surface_matrix(lat_idx, long_idx)*T*T*((4*T_arr(2) - 3*T - T_arr(3))/(2*dz)) + (3/(2*dz_arr(1)))*K_value;

                T_new = T - func/diff;
                if T_new <= 0
                    fprintf('\nNegative T in surface BC, returning NaN (lat_idx=%i, long_idx=%i, t_idx=%i)\n', lat_idx, long_idx, t_idx)
                    fprintf('%s', progress_message)
                    T_new = NaN;
                    break
                end
                difference = abs(T_new - T);
                T = T_new;
                break_iter = break_iter + 1;
                if break_iter > surface_bc_break_counter
                    fprintf('\nSurface BC iteration limit, returning NaN (lat_idx=%i, long_idx=%i, t_idx=%i)\n', lat_idx, long_idx, t_idx)
                    fprintf('%s', progress_message)
                    T_new = NaN;
                    break
                end
            end
            T_new_arr(1) = T_new;
            
            
            %% Bottom BC
            T_above = T_arr(bottom_layer_idx-1);
            T_new_arr(bottom_layer_idx) = T_above + dz_arr(bottom_layer_idx-1)*Q/(Kc_3dmat(lat_idx, long_idx, bottom_layer_idx-1)*(1 + Chi_OVER_350_POWER_3*T_above*T_above*T_above));
            
            %% Middle layers
            for layer_idx = 2:bottom_layer_idx-1
                T = T_arr(layer_idx);
                T_above = T_arr(layer_idx-1);
                alpha_var = p_arr(layer_idx)*Kc_3dmat(lat_idx, long_idx, layer_idx-1)*(1 + Chi_OVER_350_POWER_3*T_above*T_above*T_above);
                beta_var = q_arr(layer_idx)*Kc_3dmat(lat_idx, long_idx, layer_idx)*(1 + Chi_OVER_350_POWER_3*T*T*T);
                c_p_value = c0 + c1*T + c2*T*T + c3*T*T*T + c4*T*T*T*T;
                T_new_arr(layer_idx) = T + (dt/(rho_matrix(lat_idx, long_idx, layer_idx)*c_p_value))*(alpha_var*T_above - (alpha_var+beta_var)*T + beta_var*T_arr(layer_idx+1));
            end            
            
            T_new_3dmat(lat_idx, long_idx, :) = T_new_arr;
        end
    end
    T_3dmat = T_new_3dmat;
    
    if strcmp(t_wait, 'match')
        % specific ltim temperature finding routine
        T_output_3dmat(isnan(T_output_3dmat) & t_idx >= update_t_idx_matrix) = T_3dmat(isnan(T_output_3dmat) & t_idx >= update_t_idx_matrix);
    elseif t >= P*t_wait && ischar(T_mode)
        % max/min temperature finding routine
        if strcmp(T_mode, 'all')
            time_3dmat = t*ones(size(T_3dmat));
            % Extreme temperatures
            if isnan(Tmax_3dmat(1))
                Tmax_3dmat = T_3dmat;
                Tmin_3dmat = T_3dmat;
                Tmax_dtm_3dmat = time_3dmat;
                Tmin_dtm_3dmat = time_3dmat;
            else
                % max                
                cells_to_update_3dmat = T_3dmat > Tmax_3dmat;
                Tmax_3dmat = T_3dmat.*cells_to_update_3dmat + Tmax_3dmat.*~cells_to_update_3dmat;
                Tmax_dtm_3dmat = time_3dmat.*cells_to_update_3dmat + Tmax_dtm_3dmat.*~cells_to_update_3dmat;
                
                % min
                cells_to_update_3dmat = T_3dmat < Tmin_3dmat;
                Tmin_3dmat = T_3dmat.*cells_to_update_3dmat + Tmin_3dmat.*~cells_to_update_3dmat;
                Tmin_dtm_3dmat = time_3dmat.*cells_to_update_3dmat + Tmin_dtm_3dmat.*~cells_to_update_3dmat;
            end
            
            % Time series
            if t >= time_series_start && mod(t_idx, time_series_interval) == 0
                T_time_series_4dmat(:,:,:,time_series_idx) = T_3dmat;
                time_series_idx = time_series_idx + 1;
            end
            
            
        elseif strcmp(t_wait, 'match_ltim')
            switch T_mode
                case 'max'
                    cells_to_update_3dmat = T_3dmat > T_output_3dmat;
                case 'min'
                    cells_to_update_3dmat = T_3dmat < T_output_3dmat;
            end
            % matrix of hour angle for each px at previous and current time
            % step
            long_h_matrix_old = repmat(reshape(long_step_h_matrix(:,t_idx-1),1,[]), numel(lat_arr),1); 
            long_h_matrix = repmat(reshape(long_step_h_matrix(:,t_idx),1,[]), numel(lat_arr),1);
            
            % deal with 360 -> 0 wrap around
            cells_to_adjust = (long_h_matrix_old > long_h_matrix) & (long_h_matrix_old > px_measure_h_matrix);
            long_h_matrix_old(cells_to_adjust) = long_h_matrix_old(cells_to_adjust) - 360;
            cells_to_adjust = (long_h_matrix_old > long_h_matrix) & (long_h_matrix_old <= px_measure_h_matrix);
            long_h_matrix(cells_to_adjust) = long_h_matrix(cells_to_adjust) + 360;
            
            cells_to_update_3dmat = cells_to_update_3dmat...
                & (long_h_matrix_old < px_measure_h_matrix)...
                & (px_measure_h_matrix <= long_h_matrix);
            T_output_3dmat(cells_to_update_3dmat) = T_3dmat(cells_to_update_3dmat);
        else
            if isnan(T_output_3dmat(1))
                T_output_3dmat = T_3dmat;
            else
                switch T_mode
                    case 'max'
                        cells_to_update_3dmat = T_3dmat > T_output_3dmat;
                    case 'min'
                        cells_to_update_3dmat = T_3dmat < T_output_3dmat;
                end
                T_output_3dmat = T_3dmat.*cells_to_update_3dmat + T_output_3dmat.*~cells_to_update_3dmat;
            end
        end
    end
        
    if print_progress
        if mod(t_idx,print_progress) == 0 || t_idx == 1
            fraction_complete = t_idx/numel(t_arr);
            finish_dtm = datetime + (datetime - simulation_t_start_dtm)*(1-fraction_complete)/fraction_complete;
            finish_time = datestr(finish_dtm, 'HH:MM:SS dd-mmm');
            simulation_duration = datetime - simulation_t_start_dtm;
            fprintf(repmat('\b', 1, numel(progress_message)))
            progress_message = sprintf('\t\tSimulation progress: %.1f%%\n\t\tElapsed time:\t%s\n\t\tRemaining time:\t%s\n\t\tPredicted finish: %s\n',...
                100*fraction_complete, string(simulation_duration), string(finish_dtm - datetime), finish_time);
            fprintf('%s', progress_message)
        end
    end
end
% fprintf(repmat('\b', 1, numel(progress_message)))
if ~ischar(T_mode)
    T_output_3dmat = T_3dmat;
end
    
if print_progress
    simulation_duration = datetime - simulation_start_dtm;
    fprintf('\tDone in %s\n', simulation_duration)
end

if ischar(T_mode) && strcmp(T_mode, 'all')
    % Convert to proper dtm for times
    Tmax_dtm_3dmat = start_dtm + seconds(Tmax_dtm_3dmat);
    Tmin_dtm_3dmat = start_dtm + seconds(Tmin_dtm_3dmat);
    
    % Create output structure
    output = struct;
    
    % Extreme temperatures
    output.extremes.Tmax_3dmat = Tmax_3dmat;
    output.extremes.Tmax_dtm_3dmat = Tmax_dtm_3dmat;
    output.extremes.Tmin_3dmat = Tmin_3dmat;
    output.extremes.Tmin_dtm_3dmat = Tmin_dtm_3dmat;
    output.extremes.description = 'Max/min temperatures and associated datetimes for each point in the simulation. The three dimensions of the matrices give latitude, longitude and depth respectively, with the corresponding lat/long values given in metadata.crater_data and corresponding depth values given in metadata.z_arr.';
    
    % Time series
    output.time_series.T_time_series_4dmat = T_time_series_4dmat;
    output.time_series.dtm_arr = time_series_dtm_arr;
    output.time_series.time_series_interval = time_series_interval;
    output.time_series.description = 'Time series of temperatures for each point in the simulation (time resolution is time_series_interval lower than simulation time step). The four dimensions of the matrix give latitude, longitude, depth and datetime respectively, with the corresponding lat/long values given in metadata.crater_data, the corresponding depth values given in metadata.z_arr and the correspondng datetime values given in time_series_dtm_arr.';
    
    % Metadata
    output.metadata.crater_data = crater_data;
    output.metadata.z_arr = z_arr;
    output.metadata.t_wait = t_wait;
    output.metadata.scattering_power = scattering_power;
    output.metadata.num_bounces = num_bounces;
    output.metadata.custom_parameters = custom_parameters;
    output.metadata.parameters = parameters;
    output.metadata.description = 'Temperature data for location on moon, created with temperature_model_3d.';
    output.metatata.ray_tracing_datetime = ray_tracing_datetime;
    output.metadata.created = datetime;
else
    output = T_output_3dmat;
end
% save('3d_model_output_backup', 'output') % save local copy of output as backup
end