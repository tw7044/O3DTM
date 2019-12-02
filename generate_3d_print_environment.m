% GENERATE_3D_PRINT_ENVIRONMENT generates xyz coordinate matrix for 3d
% printing tsiolkovskiy crater

ppd = 16;
lat_limits = -20.40 + [-1, 1]*3.8;
long_limits = 129.1 - 0.1 + [-1, 1]*4;
[elevation_matrix, lat_arr, long_arr] = read_lola_height_data(ppd,lat_limits, long_limits);

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

elevation_matrix = elevation_matrix - min(elevation_matrix(:));

crater_data = struct;
crater_data.x_matrix = ew_matrix;
crater_data.y_matrix = ns_matrix;
crater_data.z_matrix = elevation_matrix;

save(sprintf('outputs/Tsiolkovsky_%gppd.mat', ppd), 'crater_data')

plot_3d_surface(ew_matrix, ns_matrix, elevation_matrix)
title(sprintf('ppd=%g, lat limits=[%g, %g], long limits=[%g, %g]', ppd, lat_limits(1), lat_limits(2), long_limits(1), long_limits(2)))

