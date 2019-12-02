function [T_matrix, lat_arr, long_arr, ltim_matrix, dtm_matrix] = read_diviner_temperatures(temperature, ppd)
% READ_DIVINER_TEMPERATURES loads data files containing diviner extreme
% temperatures and returns matrix of extreme temperatures. The function
% argument ('min' or 'max') defines which temperature to return.
%
% [T_matrix, lat_arr, long_arr, ltim_matrix, dtm_matrix] =
% READ_DIVINER_TEMPERATURES(temperature, ppd) returns diviner extreme
% temperatures for given temperature ('max'/'min') and ppd

if nargin == 1
    ppd = 4;
end
if strcmp(temperature, 'max')
    T_diviner = load(strcat('inputs/T_diviner/T_diviner_max_', string(ppd), 'ppd.mat'));
    T_diviner = T_diviner.T_diviner_max;
elseif strcmp(temperature, 'min')
    T_diviner = load(strcat('inputs/T_diviner/T_diviner_min_', string(ppd), 'ppd.mat'));
    T_diviner = T_diviner.T_diviner_min;
else
    error('Unexpected "temperature" value (should be "min"/"max")')
end
lat_arr = T_diviner.lat_arr;
T_matrix = T_diviner.T_diviner_matrix;
ltim_matrix = T_diviner.ltim_matrix;
long_arr = T_diviner.long_arr;
dtm_matrix = T_diviner.dtm_matrix;
end