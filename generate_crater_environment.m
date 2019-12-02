function data = generate_crater_environment(crater_name, generate_solar_angles, generate_compressed_environment, generate_compressed_solar_angles)
% GENERATE_CRATER_ENVIRONMET generates crater environment data used for 3D
% modelling. Crater names and lat/long limits are defined below.
%
% data = GENERATE_CRATER_ENVIRONMENT(...) generates crater environment,
% saves it and returns environment data file
%
% ... = GENERATE_CRATER_ENVIRONMENT(crater_name) generates_environment for
% crater_name
% 
% ... = GENERATE_CRATER_ENVIRONMENT(crater_name, generate_solar_angles)
% defines if solar angles should be calculated for crater (default=true)
% 
% ... = GENERATE_CRATER_ENVIRONMENT(crater_name, generate_solar_angles,
% generate_compressed_environment) defines if compressed environment should
% be calculated for crater (default=false)
% 
% ... = GENERATE_CRATER_ENVIRONMENT(crater_name, generate_solar_angles,
% generate_compressed_solar_angles) defines if compressed environment
% should have solar angles calculated (default=generate_solar_angles)


%% Set up
if nargin == 0
    crater_name = 'all';
end
if nargin < 2
    generate_solar_angles = true;
end
if nargin < 3
    generate_compressed_environment = false;
end
if nargin < 4
    generate_compressed_solar_angles = generate_solar_angles;
end
switch crater_name
    case 'all'
        generate_crater_environment('61n');
        generate_crater_environment('88s');
        generate_crater_environment('bruce');
        generate_crater_environment('blagg');
        generate_crater_environment('erlanger');
%         for idx = flip(1:8)
%             generate_crater_environment(sprintf('ls%d', idx));
%         end
        return
%     case 'ls'
%         for idx = flip(1:8)
%             generate_crater_environment(sprintf('ls%d', idx));
%         end
%         return
    case 'ls_4ppd'
        for idx = flip(1:8)
            generate_crater_environment(sprintf('ls%d_4ppd', idx));
        end
        return
    case 'ls_128ppd'
        for idx = flip(1:8)
            generate_crater_environment(sprintf('ls%d_128ppd', idx));
        end
        return
    case 'ls_4ppd_wide'
        for idx = flip(1:8)
            generate_crater_environment(sprintf('ls%d_4ppd_wide', idx), false, true, true);
        end
        return
    case 'ls_16ppd_wide'
        for idx = flip(1:8)
            generate_crater_environment(sprintf('ls%d_16ppd_wide', idx), false, true, true);
        end
        return
    case 'ls_128ppd_wide'
        for idx = flip(1:8)
            generate_crater_environment(sprintf('ls%d_128ppd_wide', idx), false, true, true);
        end
        return
    case 'bruce'
        lat_min = 1.0;
        lat_max = 1.4;
        long_min = 0.2;
        long_max = 0.6;
        ppd = 128;
    case 'blagg'
        lat_min = 1.09;
        lat_max = 1.34;
        long_min = 1.35;
        long_max = 1.6;
        ppd = 128;
    case 'blagg4'
        lat_min = 0.95;
        lat_max = 1.55;
        long_min = 1.20;
        long_max = 1.80;
        ppd = 4;
    case 'erlanger'
        lat_min = 86.7;
        lat_max = 87.2;
        long_min = 24.7;
        long_max = 32.7;
        ppd = 128;
    case '88s'
        lat_min = -88.45;
        lat_max = -88.14;
        long_min = 153;
        long_max = 161;
        ppd = 128;
    case '61n'
        lat_min = 61.59;
        lat_max = 61.75;
        long_min = -2.08;
        long_max = -1.80;
        ppd = 128;
%         Landing site #		latitude	longitude
%         1			-79.30		-56.00
%         2			-80.56		-37.10
%         3			-81.24		68.99
%         4			-81.35		22.80
%         5			-84.25		-4.65
%         6			-84.33		33.19
%         7			-85.33		-4.78
%         8*			-82.70		33.50
% https://quickmap.lroc.asu.edu/?layers=NrBsFYBoAZIRnpEoAsjZwLpNKG%2BscB2fDbMADlPnNAE5rDgjZWakiE2mAmOuqt3J9%2BjTOKQ9QPAohBSuiJgGZoPKgi4A6ZcuEVV6GFqLhx5OHCjckwS2iFJLJR8FXRFZJMrhqxkik5-O0NXOB9g8OVI5Tww5RcCC2UGV3VlGTTAzKSnFDYvYHV6MQk3HhyMWwzFTXNvaBTqaC0oEFV80vJVZQdSVurG6ySynuibdrVU-rMG0Jnuxs9jWbdGypX6t2lGatBllp5V2PGF71Rm1u6IS7NRxsEZwdjbxYTX0aI%2BqvaiU9YtFZFjxvsYgQ0jpdwWspFDjmpEgCgZ84rJfsNCglHpiiNNMRQCuw3BRlioKBsyf8icoKKCsGV8oSEO1RK5VDwqSo1BStrFEXJtqjMaBsUT1CxgrE8dSJWygnKeQyQbsQCgOZc9Eg1RjjJrgGqhcY6FsUOBObZTTyteAdS07os4NLmWtoIaudA6YtGpLoL7guAeEz9KBOcGdbwKBcsrhgsVRfTyOArLJnUnPUgk-ymEn44nHf6eAcTImKkYWkRgyRNNTfVaXeHzJggA&probeTool=33.51195547%2C-82.70705462&extent=-120.98977689607486%2C-53.993956418505554%2C120.73450462459725%2C-50.63083898119902&proj=17
    % *** ADD NEW CRATER ENVIRONMENTS HERE ***
    otherwise
        if strcmp(crater_name(1:2), 'ls')
            %%%need to add if double or single digit?
            %%is crater_name(4) a digit?
            if(~isnan(str2double('crater_name(4)')) == 1)
                crater_numel = str2double(carter_name(3:4));
            else
                crater_numel = str2double(carter_name(3));
            end
            switch crater_numel
                case '1'
                    ls_lat = -79.30;
                    ls_long = -56.00;
                case '2'
                    ls_lat = -80.56;
                    ls_long = -37.10;
                case '3'
                    ls_lat = -81.24;
                    ls_long = 68.99;
                case '4'
                    ls_lat = -81.35;
                    ls_long = 22.80;
                case '5'
                    ls_lat = -84.25;
                    ls_long = -4.65;
                case '6'
                    ls_lat = -84.33;
                    ls_long = 33.19;
                case '7'
                    ls_lat = -85.33;
                    ls_long = -4.78;
                case '8'
                    ls_lat = -82.70;
                    ls_long = 33.50;
                otherwise
                    error('Unknown crater')
            end
            if strcmp(crater_name(end-4:end), '_4ppd')
               lat_buffer = 4;
               long_buffer = 4;
               ppd = 4;
            elseif strcmp(crater_name(end-6:end), '_128ppd')
                lat_buffer = 0.15;
                long_buffer = 0.15;
                ppd = 128;
            elseif strcmp(crater_name(end-11:end), '_128ppd_wide')
                lat_buffer = 0.15;
                long_buffer = 0.75;
                ppd = 128;
            elseif strcmp(crater_name(end-10:end), '_16ppd_wide')
                lat_buffer = 2;
                long_buffer = 10;
                ppd = 16;
            elseif strcmp(crater_name(end-9:end), '_4ppd_wide')
                lat_buffer = 4;
                long_buffer = 20;
                ppd = 4;
            end
            lat_min = ls_lat - lat_buffer;
            lat_max = ls_lat + lat_buffer;
            long_min = ls_long - long_buffer;
            long_max = ls_long + long_buffer;
        else
            error('Unknown crater')
        end
end
%% Load LOLA data
fprintf('Loading LOLA elevation data... ')
[elevation_matrix, lat_arr, long_arr] = read_lola_height_data(ppd, [lat_min, lat_max], [long_min, long_max]);

% get distances
ref_sphere = referenceSphere('moon');
center_lat = mean(lat_arr);
center_long = mean(long_arr);
ns_matrix = NaN(size(elevation_matrix));
ew_matrix = NaN(size(elevation_matrix));
for lat_idx = 1:numel(lat_arr)
    lat = lat_arr(lat_idx);
    for long_idx = 1:numel(long_arr)
        long = long_arr(long_idx);
        [xEast,yNorth,~] = geodetic2enu(lat,long,0,center_lat,center_long,0,ref_sphere);
        ew_matrix(lat_idx, long_idx) = xEast;
        ns_matrix(lat_idx, long_idx) = yNorth;
    end
end

fprintf('Done\n')

%% Get A0 data
fprintf('Calculating albedo data... ')
if ppd > 4
    albedo_ppd = 8;
else
    albedo_ppd = 4;
end
[albedo_matrix, albedo_lat_arr, albedo_long_arr] = read_lola_A0_data(albedo_ppd); % use highest resolution albedo data possible
parameters = define_parameters;
albedo_matrix = scale_data(albedo_matrix, parameters.A0);
[albedo_long_matrix, albedo_lat_matrix] = meshgrid(albedo_long_arr, albedo_lat_arr);
[long_matrix, lat_matrix] = meshgrid(long_arr, lat_arr);
albedo_matrix = interp2(albedo_long_matrix, albedo_lat_matrix, albedo_matrix, long_matrix, lat_matrix);
A0 = mean(albedo_matrix(:));
fprintf('Done\n')

%% Load Diviner data
fprintf('Loading Diviner temperature data...\n')
[Tmax_matrix, Tmin_matrix, Tmax_ltim_matrix, Tmin_ltim_matrix, Tmax_dtm_matrix, Tmin_dtm_matrix] = find_diviner_extreme_temperatures(ppd, lat_arr, long_arr);


%% Save data
fprintf('Generating environment for "%s"... ', crater_name)
data = struct;
data.crater_name = crater_name;
data.lat_arr = lat_arr;
data.long_arr = long_arr;
data.ew_matrix = ew_matrix;
data.ns_matrix = ns_matrix;
data.center_lat = mean(lat_arr);
data.center_long = mean(long_arr);
data.elevation_matrix = elevation_matrix;
data.Tmax_matrix = Tmax_matrix;
data.Tmax_ltim_matrix = Tmax_ltim_matrix;
data.Tmax_dtm_matrix = Tmax_dtm_matrix;
data.Tmin_matrix = Tmin_matrix;
data.Tmin_ltim_matrix = Tmin_ltim_matrix;
data.Tmin_dtm_matrix = Tmin_dtm_matrix;
data.A0 = A0;
data.ppd = ppd;
data.description = 'Data for lunar crater';
data.created = datetime;

target_name = create_static_path(sprintf('crater_environments/%s.mat', crater_name));
save(target_name, 'data')
fprintf('Done\n')
if generate_solar_angles
    generate_crater_ray_tracing(crater_name, false);
    generate_crater_ray_tracing(crater_name, true);
end
if generate_compressed_environment
    result = generate_compressed_crater_environment(crater_name);
    if isfield(result, 'compression_ratio') && (generate_solar_angles || generate_compressed_solar_angles)
        crater_name = sprintf('%s_compressed', crater_name);
        generate_crater_ray_tracing(crater_name, false);
        generate_crater_ray_tracing(crater_name, true);
    end
end
end