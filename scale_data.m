function [output_matrix, scaling_factor] = scale_data(raw_data_matrix, desired_mean, scaling_factor)
% SCALE_DATA scales raw_data_matrix so that the mean of the output_matrix
% is desired_mean.
% 
% [output_matrix, scaling_factor] = SCALE_DATA(...) returns scaled data and
% scaling factor
%
% SCALE_DATA(raw_data_matrix, desired_mean) scales raw_data_matrix to
% desired_mean
%
% SCALE_DATA(raw_data_matrix, desired_mean, scaling_factor) scales
% raw_data_matrix by scaling factor

if nargin < 3
    raw_mean = nanmean(nanmean(raw_data_matrix));
    scaling_factor = desired_mean/raw_mean;
end
output_matrix = raw_data_matrix.*scaling_factor;
end