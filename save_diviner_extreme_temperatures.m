function save_diviner_extreme_temperatures(ppd)
% SAVE_DIVINER_EXTREME_TEMPERATURES saves data files for maximum and
% minimum diviner temperatures for given ppd
%
% SAVE_DIVINER_EXTREME_TEMPERATURES(ppd) specifies the resolution to save
% extreme temperatures for

if nargin == 0
    ppd = 4;
end
[Tmax_matrix, Tmin_matrix, Tmax_ltim_matrix, Tmin_ltim_matrix, Tmax_dtm_matrix, Tmin_dtm_matrix, lat_arr, long_arr] = ...
    find_diviner_extreme_temperatures(ppd);

fprintf('Saving data...\n')

T_diviner_max = struct;
T_diviner_max.T_diviner_matrix = Tmax_matrix;
T_diviner_max.ltim_matrix = Tmax_ltim_matrix;
T_diviner_max.dtm_matrix = Tmax_dtm_matrix;
T_diviner_max.lat_arr = lat_arr;
T_diviner_max.long_arr = long_arr;
T_diviner_max.created = datetime;
T_diviner_max.ppd = ppd;
T_diviner_max.description = 'Maximum Diviner temperature data for each pixel, calculated using find_diviner_extreme_temperatures.m';
target_name = strcat('inputs/T_diviner/T_diviner_max_', string(ppd), 'ppd.mat');
save(target_name, 'T_diviner_max')

T_diviner_min = struct;
T_diviner_min.T_diviner_matrix = Tmin_matrix;
T_diviner_min.ltim_matrix = Tmin_ltim_matrix;
T_diviner_min.dtm_matrix = Tmin_dtm_matrix;
T_diviner_min.lat_arr = lat_arr;
T_diviner_min.long_arr = long_arr;
T_diviner_min.created = datetime;
T_diviner_min.ppd = ppd;
T_diviner_min.description = 'Minimum Diviner temperature data for each pixel, calculated using find_diviner_extreme_temperatures.m';
target_name = strcat('inputs/T_diviner/T_diviner_min_', string(ppd), 'ppd.mat');
save(target_name, 'T_diviner_min')

fprintf('Done\n')
end