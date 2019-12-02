function [dtm_arr, sub_sun_long_arr, decl_arr] = read_horizons_data
% READ_HORIZONS_DATA loads data file for JPL horizons data (defining lunar
% rotation and seasonal variations).
%
% [dtm_arr, sub_sun_long_arr, decl_arr] = READ_HORIZONS_DATA

horizons_data = load(create_static_path('ephemerides/horizons_data.mat'));
horizons_data = horizons_data.horizons_data;
dtm_arr = horizons_data.dtm_arr;
sub_sun_long_arr =  horizons_data.sub_sun_long_arr;
decl_arr = horizons_data.decl_arr;
end