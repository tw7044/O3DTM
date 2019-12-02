function data = generate_simulated_crater_environment(center_lat, generate_solar_angles, n_pts, crater_depth, crater_radius, crater_power)
% GENERATE_SIMULATED_CRATER_ENVIRONMENT generates a simulated bowl shaped
% crater environment with given parameters and saves data file for crater
% 
% data = GENERATE_SIMULATED_CRATER_ENVIRONMENT(...) returns crater data
% file
%
% GENERATE_SIMULATED_CRATER_ENVIRONMENT(center_lat) defines central
% latitude of crater
%
% GENERATE_SIMULATED_CRATER_ENVIRONMENT(center_lat, generate_solar_angles)
% controls if ray tracing should be carried out for crater
%
% GENERATE_SIMULATED_CRATER_ENVIRONMENT(center_lat, generate_solar_angles,
% n_pts) defines number of points in simulation grid
%
% GENERATE_SIMULATED_CRATER_ENVIRONMENT(center_lat, generate_solar_angles,
% n_pts, crater_depth) defines depth of crater 
%
% GENERATE_SIMULATED_CRATER_ENVIRONMENT(center_lat, generate_solar_angles,
% n_pts, crater_depth, crater_radius) defines radius of crater
%
% GENERATE_SIMULATED_CRATER_ENVIRONMENT(center_lat, generate_solar_angles,
% n_pts, crater_depth, crater_radius, crater_power) defines steepness of
% crater walls

if nargin == 0
    generate_simulated_crater_environment(0);
    generate_simulated_crater_environment(30);
    generate_simulated_crater_environment(60);
    generate_simulated_crater_environment(85);
    generate_crater_ray_tracing('seasons')
    return
end
if nargin < 2
    generate_solar_angles = false;
end
if nargin < 3
    n_pts = 16;
end
if nargin < 4
    crater_depth = 500;
end
if nargin < 5
    crater_radius = 1e3;
end
if nargin < 6
    crater_power = 2;
end
ref_sphere = referenceSphere('moon');
center_long = 0;
center_lat = abs(center_lat);

map_radius = ((n_pts + 2)/n_pts)*crater_radius;

range_arr = linspace(-map_radius, map_radius, n_pts+2);
lat_arr = aer2geodetic(0, 0, range_arr, center_lat, center_long, 0, ref_sphere);
[~, long_arr] = aer2geodetic(90, 0, range_arr, center_lat, center_long, 0, ref_sphere);

elevation_matrix = NaN(n_pts+2);
ns_matrix = NaN(size(elevation_matrix));
ew_matrix = NaN(size(elevation_matrix));
for lat_idx = 1:numel(lat_arr)
    lat = lat_arr(lat_idx);
    for long_idx = 1:numel(long_arr)
        long = long_arr(long_idx);
        [xEast,yNorth,~] = geodetic2enu(lat,long,0,center_lat,center_long,0,ref_sphere);
        ew_matrix(lat_idx, long_idx) = xEast;
        ns_matrix(lat_idx, long_idx) = yNorth;
        
        if norm([xEast,yNorth]) > crater_radius
            elevation_matrix(lat_idx, long_idx) = crater_depth;
        else
            elevation = crater_depth*(norm([xEast,yNorth])/crater_radius)^crater_power;
            elevation_matrix(lat_idx, long_idx) = elevation;
        end
    end
end
ref_sphere = referenceSphere('moon');
elevation_matrix = elevation_matrix + ref_sphere.Radius;

crater_name = sprintf('simulated_%gn', center_lat);
data = struct;
data.crater_name = crater_name;
data.lat_arr = lat_arr;
data.long_arr = long_arr;
data.ew_matrix = ew_matrix;
data.ns_matrix = ns_matrix;
data.center_lat = center_lat;
data.center_long = center_long;
data.elevation_matrix = elevation_matrix;
data.simulated_crater = true;
data.crater_depth = crater_depth;
data.crater_radius = crater_radius;
data.crater_power = crater_power;
data.description = 'Data for simulated lunar crater' ;
data.created = datetime;

target_name = create_static_path(sprintf('crater_environments/%s.mat', crater_name));
save(target_name, 'data')

if generate_solar_angles
    generate_crater_ray_tracing(crater_name, false);
    generate_crater_ray_tracing(crater_name, true);
end
end