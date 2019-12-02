function data = plot_compare_model_diviner_temperatures(data_path, plot_TT, color, keep_old_plot)
% PLOT_COMPARE_MODEL_DIVINER_TEMPERATURES plots comparison of 1D model and
% Diviner temperatures for identifying where 1D model fails.
%
% data = PLOT_COMPARE_MODEL_DIVINER_TEMPERATURES(...) returns data file of
% model temperatures
%
% PLOT_COMPARE_MODEL_DIVINER_TEMPERATURES(data_path) defines path to data
% file containing model extreme temperatures
% 
% PLOT_COMPARE_MODEL_DIVINER_TEMPERATURES(data_path, plot_TT) switches
% between two different plot types (plot_TT = true/false)
%
% PLOT_COMPARE_MODEL_DIVINER_TEMPERATURES(data_path, plot_TT, color)
% defines colors to use for plotting
%
% PLOT_COMPARE_MODEL_DIVINER_TEMPERATURES(data_path, plot_TT, color,
% keep_old_plot) defines if old plot should be cleared

if nargin < 2
    plot_TT = true;
end
if nargin < 3
    keep_old_plot = false;
end
data = load(data_path);
data = data.data;
Tdiviner_matrix = read_diviner_temperatures(data.T_mode);
[~, lat_matrix] = meshgrid(data.long_arr, data.lat_arr);
if ~keep_old_plot
    clf
end
hold on
box on

if plot_TT
    marker_size = 1;
    color_matrix = abs(lat_matrix);
    color_matrix = color_matrix(:);
    if nargin == 3
        color_matrix = color;
    end
    colormap(flipud(viridis))
    xlabel('Diviner')
    ylabel('Model')
    if nargin < 3
        caxis([0,90])
        c = colorbar;
        c.Label.String = sprintf('Latitude (%sN / %sS)', char(176), char(176));
    end
    scatter(Tdiviner_matrix(:), data.Tsimulated_matrix(:), marker_size, color_matrix, '.', 'MarkerEdgeAlpha',0.3)
    lim = [min(min(Tdiviner_matrix(:)), min(data.Tsimulated_matrix(:))), max(max(Tdiviner_matrix(:)), max(data.Tsimulated_matrix(:)))];
    plot(lim, lim, 'k');%'Color', [0.75, 0.75, 0.75])
    xlim(lim)
    ylim(lim)
else
    marker_size = 5;
    x_arr = abs(lat_matrix(:));
    y_arr = (data.Tsimulated_matrix(:) - Tdiviner_matrix(:));%./Tdiviner_matrix(:);
    c_arr = Tdiviner_matrix(:);
    x_bins = min(x_arr):0.25:max(x_arr);
    y_bins = (round(min(y_arr))-0.5):0.5:(round(max(y_arr))+0.5);
    if nargin >= 3
        c_arr = color;
    else
        c = colorbar;
        c.Label.String = 'T_{Diviner} (K)'; 
    end
    xlim([0,90])
    xticks(0:15:90)
%     scatter(x_arr, y_arr, marker_size, c_arr, '.', 'MarkerEdgeAlpha',0.05)
    histogram2(x_arr, y_arr, x_bins, y_bins,'DisplayStyle','tile');
    xlabel(sprintf('Latitude (%sN / %sS)', char(176), char(176)))
    ylabel('T_{Model} - T_{Diviner} (K)')
    plot([0,90], [0,0], 'k--');
end

end