function varargout = plot_diff_map(lat_arr, long_arr, value_matrix, title_str, cbar_str)
% PLOT_DIFF_MAP plots map of difference between two different maps
% (i.e. version of plot_map for comparing different parameter maps)
% Use two arguments for the filepaths to the two different data files (from
% generate_map) to compare, or use 3+ arguments for plotting a specific
% matrix.
% 
% varargout = PLOT_DIFF_MAP(lat_arr, long_arr, value_matrix, title_str, cbar_str)

if nargin < 4
    title_str = '';
end
if nargin < 5
    cbar_str = '';
end
if nargin == 2
    data_path = lat_arr;
    data_from_file = load(data_path);
    data = data_from_file.data;
    value_matrix_1 = data.value_matrix;
    title_str_1 = data.plot_title;
    
    data_path = long_arr;
    data_from_file = load(data_path);
    data = data_from_file.data;
    value_matrix_2 = data.value_matrix;
    title_str_2 = data.plot_title;
    
    lat_arr = data.lat_arr;
    long_arr = data.long_arr;
    cbar_str = sprintf('%s difference', data.plot_label);

    value_matrix = value_matrix_1 - value_matrix_2;
    title_str = sprintf('%s - %s', title_str_1{1}, title_str_2{1});
end


plot_h = plot_map(lat_arr, long_arr, value_matrix, title_str, cbar_str);
divergent_colormap

if nargout > 0
    [varargout{1:nargout}] = plot_h;
end

end