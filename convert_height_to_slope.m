function [aspect_matrix, slope_matrix, lat_matrix, long_matrix] = convert_height_to_slope(radius_matrix, lat_arr, long_arr)
% CONVERT_HEIGHT_TO_SLOPE converts lunar radius matrix to aspect and slope
% matrices.
% 
% [aspect_matrix, slope_matrix, lat_matrix, long_matrix] =
% CONVERT_HEIGHT_TO_SLOPE(radius_matrix, lat_arr, long_arr) returns the
% slope aspect matrix, slope magnitude matrix, latitude matrix and
% longitude matrix for given matrix of lunar radii and corresponding
% lat/long arrays.

lat_matrix = repmat(reshape(lat_arr,[],1), 1, numel(long_arr));
long_matrix = repmat(long_arr, numel(lat_arr), 1);
ref_sphere = referenceSphere('moon');
radius_matrix = radius_matrix - ref_sphere.Radius;
[aspect_matrix, slope_matrix] = gradientm(lat_matrix, long_matrix, radius_matrix, ref_sphere);
aspect_matrix = aspect_matrix - 180; % redefine aspect to version needed for theta formula
aspect_matrix(isnan(aspect_matrix)) = 0; % ensure doesn't break for slope = 0
end