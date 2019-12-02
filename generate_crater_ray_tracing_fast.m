function data = generate_crater_ray_tracing_fast(crater_name, use_seasons)
% GENERATE_CRATER_RAY_TRACING_FAST test version of generate_crater_ray
% tracing to test performance improvements

warning('This is a test function, use generate_crater_ray_testing for real version')

if nargin == 0 || strcmp(crater_name, 'all')
    generate_crater_ray_tracing_fast('bruce', false);
    generate_crater_ray_tracing_fast('blagg', false);
    generate_crater_ray_tracing_fast('erlanger_compressed', false);
    generate_crater_ray_tracing_fast('bruce', true);
    generate_crater_ray_tracing_fast('blagg', true);
    generate_crater_ray_tracing_fast('erlanger_compressed', true);
    return
end
if strcmp(crater_name, 'seasons')
    crater_name_arr = {'simulated_0n', 'simulated_85n', 'simulated_60n', 'simulated_30n'};
    for crater_idx = 1:numel(crater_name_arr)
        crater_name = crater_name_arr{crater_idx};
        generate_crater_ray_tracing_fast(crater_name);
    end
    for crater_idx = 1:numel(crater_name_arr)
        crater_name = crater_name_arr{crater_idx};
        generate_crater_ray_tracing_fast(crater_name, true);
    end
end
if strcmp(crater_name(1:3), 'ls_')
    for ls_idx = flip(1:8)
        generate_crater_ray_tracing_fast(sprintf('ls%d%s', ls_idx, crater_name(3:end)), false);
        generate_crater_ray_tracing_fast(sprintf('ls%d%s', ls_idx, crater_name(3:end)), true);        
    end
    return
end
if nargin < 2
    generate_crater_ray_tracing_fast(crater_name, false);
    generate_crater_ray_tracing_fast(crater_name, true);
end

if use_seasons
    seasons_str = 'seasonal';
else
    seasons_str = 'non-seasonal';
end
tic
fprintf('Calculating %s pixel information for "%s"\n   ', seasons_str, crater_name)

%% Process data
source_data = load_crater_environment(crater_name);

lat_arr = source_data.lat_arr;
long_arr = source_data.long_arr;
elevation_matrix = source_data.elevation_matrix;

parameters = define_parameters;
P = parameters.P;
const_decl = parameters.decl;
dt = parameters.dt;
phase_start_dtm = parameters.phase_start_dtm;
phase_end_dtm = parameters.phase_end_dtm;

%% Get info from inputs
% Get information about terrain and linear distances etc.
ref_sphere = referenceSphere('moon');
r_moon = ref_sphere.Radius; % moon radius in m
elevation_matrix = elevation_matrix - r_moon; % make elevations relative to lunar radius
raster_ref = georefcells([min(lat_arr), max(lat_arr)], [min(long_arr), max(long_arr)], size(elevation_matrix), 'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');

[aspect_matrix, slope_matrix] = convert_height_to_slope(elevation_matrix, lat_arr, long_arr);
[~,~,max_distance] = geodetic2aer(...
    lat_arr(1), long_arr(1), min(reshape(elevation_matrix,[],1)),...
    lat_arr(end), long_arr(end), max(reshape(elevation_matrix,[],1)),...
    ref_sphere);


%% Generate larger grid for calculating shadowing
% Slightly larger grid ensures floating point errors do not break shadowing
% calculations at the edge of the grid
extended_elevation_matrix = zeros(size(elevation_matrix) + [2,2]);
extended_elevation_matrix(2:end-1,2:end-1) = elevation_matrix;
extended_elevation_matrix(2:end-1,1) = elevation_matrix(:,1);
extended_elevation_matrix(2:end-1,end) = elevation_matrix(:,end);
extended_elevation_matrix(1,2:end-1) = elevation_matrix(1,:);
extended_elevation_matrix(end,2:end-1) = elevation_matrix(end,:);
extended_elevation_matrix(1,1) = elevation_matrix(1,1);
extended_elevation_matrix(1,end) = elevation_matrix(1,end);
extended_elevation_matrix(end,1) = elevation_matrix(end,1);
extended_elevation_matrix(end,end) = elevation_matrix(end,end);

d_lat = abs(lat_arr(1) - lat_arr(2));
d_long = abs(long_arr(1) - long_arr(2));
extended_raster_ref = georefcells([min(lat_arr)-d_lat, max(lat_arr)+d_lat], [min(long_arr)-d_long, max(long_arr)+d_long], size(extended_elevation_matrix), 'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');

%% Prepare declination data
% If using seasons, process declination data for use in simulation
if use_seasons
    long = mean(long_arr);
    [raw_dtm_arr, sub_sun_long_arr, decl_arr] = read_horizons_data();
    long_interp_arr = sub_sun_long_arr((raw_dtm_arr > phase_start_dtm) & (raw_dtm_arr < phase_end_dtm));
    dtm_interp_arr = raw_dtm_arr((raw_dtm_arr > phase_start_dtm) & (raw_dtm_arr < phase_end_dtm));
    % ensure both have consistent ranges
    long = mod(long + 180, 360) - 180;
    long_interp_arr = mod(long_interp_arr + 180, 360) - 180;
    [long_interp_arr, unique_idx_arr] = unique(long_interp_arr);
    dtm_interp_arr = dtm_interp_arr(unique_idx_arr);
    start_dtm = interp1(long_interp_arr(2:end-1), dtm_interp_arr(2:end-1), long); % start simulation at local midday
    end_dtm = raw_dtm_arr(end);
    dtm_arr = start_dtm:seconds(dt):end_dtm;
    decl_arr = interp1(raw_dtm_arr, decl_arr, dtm_arr);
    sub_sun_long_arr = interp1(raw_dtm_arr, sub_sun_long_arr, dtm_arr);
    h_arr = long - sub_sun_long_arr;
    h_arr = mod(h_arr + 180, 360) - 180;
else
    h_arr = [];
    decl_arr = [];
end
decl_arr = single(decl_arr);
h_arr = single(h_arr);
clear raw_dtm_arr dtm_interp_arr sub_sun_long_arr

%% Initialise variables
if use_seasons
    t_steps = numel(decl_arr);
else
    t_steps = P/dt;
end
dh = 360/t_steps;
theta_3dmat = NaN(numel(lat_arr), numel(long_arr), t_steps, 'single');
vis_4dmat = false(numel(lat_arr), numel(long_arr), numel(lat_arr), numel(long_arr));

%% Calculate values

progress_idx = 0;
calculation_t_start_dtm = datetime;
progress_message = '';
for lat_idx = 1:numel(lat_arr)
    lat = lat_arr(lat_idx);

    %% Prepare for parfor iteration
    % Set up local copies so parfor works
    for long_idx = 1:numel(long_arr)
        long = long_arr(long_idx);

        %% Get information about current location
        slope1 = slope_matrix(lat_idx, long_idx);
        aspect1 = aspect_matrix(lat_idx, long_idx);
        height_px = elevation_matrix(lat_idx, long_idx);
        if lat >= 0
            aspect = aspect1;
            slope = slope1;
        else
            aspect = -aspect1;
            slope = -slope1;
        end
        
        if use_seasons
            h_arr_local = h_arr + long - mean(long_arr);
            decl_arr_local = decl_arr;
        else
            h_arr_local = [];
            decl_arr_local = [];
        end
        
        %% Generate visibility matrix
        vis_4dmat(lat_idx, long_idx, :, :) = viewshed(elevation_matrix, raster_ref, lat, long, 0, 0, 'AGL', 'AGL', r_moon);

        %% Generate solar incidence angles & solar fluxes
        % Generate solar incidence angles using formula from Braun and Mitchell [1983]
        t_calculation_step_length = 1e3; % calculate each 1e3 time steps at a time
        t_calculation_steps = ceil(t_steps/t_calculation_step_length);
        theta_matrix = NaN(t_calculation_steps, t_calculation_step_length, 'single');
        parfor t_calculation_step_idx = 1:t_calculation_steps
            t_arr_start = 1+(t_calculation_step_idx-1)*t_calculation_step_length;
            t_arr_end = t_arr_start + t_calculation_step_length-1;
            if t_arr_end > t_steps
                t_arr_end = t_steps;
            end
            t_arr = t_arr_start:t_arr_end;
            if numel(t_arr) == 0
                continue
            end
            
            if use_seasons
                decl = decl_arr_local(t_arr);
                h = h_arr_local(t_arr);
            else
                decl = const_decl*ones(size(t_arr));
                h = dh*t_arr - dh;
            end
            h = mod(h + 180, 360) - 180;
            sigma_ew = ones(size(t_arr), 'single');
            sigma_ew(abs(h) > acosd(cotd(lat).*tand(decl))) = -1;
            sigma_ns = ones(size(t_arr), 'single');
            sigma_ns(lat*(lat-decl) < 0) = -1;
            sigma_w = ones(size(t_arr), 'single');
            sigma_w(h < 0) = -1;
            theta_z = acos(sind(decl).*sind(lat) + cosd(decl).*cosd(lat).*cosd(h));
            gamma_so = real(asind(sind(h).*cosd(decl)./sin(theta_z)));
            gamma_so(theta_z == 0) = 0; % Avoid division by 0 error
            % gamma_so may end up complex if it is ~90deg as rounding errors in
            % the asind(...) function can make its argument slightly larger
            % than 1, making gamma_so slightly complex. Therefore, only use the
            % real part of gamma_so to avoid the tiny complex part affecting
            % later calculations.
            gamma_s = sigma_ew.*sigma_ns.*gamma_so + ((1-sigma_ew.*sigma_ns)/2).*sigma_w*180;
            theta = acos(cos(theta_z).*cosd(slope) + sin(theta_z).*sind(slope).*cosd(gamma_s - aspect));
            theta(theta_z > pi/2) = NaN;
            theta(theta > pi/2) = NaN;
            
            az_sun = NaN(size(t_arr));
            az_sun(lat >= decl) = 180+gamma_s(lat >= decl);
            az_sun(lat < decl) = -gamma_s(lat < decl);
            elev_sun = 90-rad2deg(theta_z);
            slant_range = 1.1*max_distance; % ensure outside of grid
            [lat_sun,long_sun,height_sun] = aer2geodetic(double(az_sun),double(elev_sun),slant_range,lat,long,height_px,ref_sphere);
            
            sun_vis = los2(extended_elevation_matrix, extended_raster_ref, lat*ones(size(lat_sun)), long*ones(size(lat_sun)), lat_sun, long_sun, 0, height_sun, 'AGL', 'MSL', r_moon);
            sun_vis(theta_z == 0) = true;
            theta(~sun_vis) = NaN;
            if numel(theta) < t_calculation_step_length
                theta(end+1:t_calculation_step_length) = NaN;
            end
            theta_matrix(t_calculation_step_idx, :) = theta;
        end
        theta_arr = [];
        for idx = 1:t_calculation_steps
            % get values in correct order for time series
            theta_arr(end+1:end+t_calculation_step_length) = theta_matrix(idx,:);
        end
        theta_3dmat(lat_idx,long_idx,:) = theta_arr(1:t_steps);
        
        %% Print progress
        progress_idx = progress_idx + 1;
        fraction_complete = progress_idx/(numel(lat_arr)*numel(long_arr));
        finish_dtm = datetime + (datetime - calculation_t_start_dtm)*(1-fraction_complete)/fraction_complete;
        finish_time = datestr(finish_dtm, 'HH:MM:SS dd-mmm');
        simulation_duration = datetime - calculation_t_start_dtm;
        fprintf(repmat('\b', 1, numel(progress_message)))
        progress_message = sprintf('\tCalculation progress: %.1f%%\n\tElapsed time:\t%s\n\tRemaining time:\t%s\n\tPredicted finish: %s\n',...
            100*fraction_complete, string(simulation_duration), string(finish_dtm - datetime), finish_time);
        fprintf('%s', progress_message)
    end
end

%% Save data
data = struct;
data.crater_name = source_data.crater_name;
data.lat_arr = source_data.lat_arr;
data.long_arr = source_data.long_arr;
data.ew_matrix = source_data.ew_matrix;
data.ns_matrix = source_data.ns_matrix;
data.center_lat = source_data.center_lat;
data.center_long = source_data.center_long;
data.elevation_matrix = source_data.elevation_matrix;
data.use_seasons = use_seasons;
data.theta_3dmat = theta_3dmat;
data.dt = dt;
data.vis_4dmat = vis_4dmat;
if use_seasons
    data.dtm_arr = dtm_arr;
    data.h_arr = h_arr;
    data.description = 'Solar angle and pixel visibility data for lunar crater at times given by dtm_arr';
    file_season_str = 'seasons';
else
    data.description = 'Solar angle and pixel visibility data for lunar crater through lunar day';
    file_season_str = 'no_seasons';
end
data.created = datetime;
target_name = create_static_path(sprintf('crater_environments/ray_tracing/%s/%s_ray_tracing.mat', file_season_str, crater_name));
save(target_name, 'data', '-v7.3') % save as v7.3 due to big file size
fprintf('Done\n')
toc
end