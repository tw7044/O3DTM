% ANALYSE_LS_DIVINER_TEMPERATURES containes various routines to alalyse
% landing site diviner temperatures

ls_name = 'ls8_4ppd_wide_compressed';
output = load(create_static_path(sprintf('crater_environments/diviner_temperatures/%s_diviner_temperatures.mat', ls_name)));
output = output.output;

ls_lat_arr = output.crater_data.lat_arr;
ls_long_arr = output.crater_data.long_arr;
% [long_matrix, lat_matrix] = meshgrid(ls_long_arr, ls_lat_arr);


data_lat_arr = output.lat_arr;
data_long_arr = output.long_arr;

if ~exist('data', 'var') || ~isstruct(data) || ~isfield(data, 'metadata') || ~isfield(data.metadata, 'crater_data') || ~isfield(data.metadata.crater_data, 'crater_name') || ~strcmp(data.metadata.crater_data.crater_name, output.crater_data.crater_name) || numel(data.metadata.crater_data.long_arr) ~= numel(output.crater_data.long_arr)
    fprintf('Loading data...\n')
    load(create_static_path(sprintf('outputs/ls_temperatures/%s_temperatures.mat', ls_name)));
    if data.metadata.created < datetime(2018, 8, 28, 0, 0, 0)
        warning('Old data file from before albedo fix (created %s)', data.metadata.created)
    end
end

ls_name = strrep(ls_name, '_', '\_');

%% Choose mode
diviner_consistency = 0;
surface_map = 0;
single_point = 1;
layer_time_series = 0;
layer_depth_profile = 0;
profile_3d_old = 0;
profile_3d = 0;
stable_condition_T = 0;
stable_condition_z = 0;
plot_ls_locations = 0;
compare_deep_surface_T = 0;


%% Diviner consistency
if diviner_consistency
    plot_matrix = NaN(size(output.crater_data.elevation_matrix));
    for lat_idx = 1:numel(output.crater_data.lat_arr)
        fprintf('.')
        lat = output.crater_data.lat_arr(lat_idx);
        for long_idx = 1:numel(output.crater_data.long_arr)
            long = output.crater_data.long_arr(long_idx);
            bool_arr = output.lat_arr == lat & output.long_arr == long & output.dtm_arr >= min(data.time_series.dtm_arr) & output.dtm_arr <= max(data.time_series.dtm_arr);
            T_arr = output.T_arr(bool_arr);
            dtm_arr = output.dtm_arr(bool_arr);
            err_value = 0;
            T_model_arr = squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, 1, :));
            for T_idx = 1:numel(T_arr)
                T_model = interp1(data.time_series.dtm_arr, T_model_arr, dtm_arr(T_idx));
                err_value = err_value + (T_model - T_arr(T_idx))^2;
            end
            plot_matrix(lat_idx, long_idx) = (err_value/numel(T_arr))^0.5;
        end
    end
    plot_3d_surface(output.crater_data.ew_matrix, output.crater_data.ns_matrix, output.crater_data.elevation_matrix, plot_matrix, ls_name, 'Model temperature error')
    shading faceted
end


%% Surface map
if surface_map
    T_arr = output.T_arr;
%     number_matrix = NaN(numel(ls_lat_arr), numel(ls_long_arr));
%     for lat_idx = 1:numel(ls_lat_arr)
%         lat = ls_lat_arr(lat_idx);
%         for long_idx = 1:numel(ls_long_arr)
%             long = ls_long_arr(long_idx);
%             number_matrix(lat_idx,long_idx) = numel(T_arr(data_lat_arr == lat & data_long_arr == long));
%         end
%     end
    plot_matrix = data.extremes.Tmax_3dmat(:,:,1) - output.crater_data.Tmax_matrix;
    plot_3d_surface(output.crater_data.ew_matrix, output.crater_data.ns_matrix, output.crater_data.elevation_matrix, plot_matrix, ls_name);% 'Number of Diviner measurements per pixel')
    shading faceted
    divergent_colormap
end


%% Single point
if single_point
%%%%ls8: (-82.70, 33.50)
    lat_idx = round(numel(output.crater_data.lat_arr)/2);
    long_idx = round(numel(output.crater_data.long_arr)/2);
    
    lat = output.crater_data.lat_arr(lat_idx);
    long = output.crater_data.long_arr(long_idx);
    
    T_diviner_arr = output.T_arr(data_lat_arr == lat & data_long_arr == long);
    ltim_arr = output.ltim_arr(data_lat_arr == lat & data_long_arr == long);
    jd_arr = output.jd_arr(data_lat_arr == lat & data_long_arr == long);
    start_dtm_arr = output.start_dtm_arr(data_lat_arr == lat & data_long_arr == long);
    end_dtm_arr = output.end_dtm_arr(data_lat_arr == lat & data_long_arr == long);
    
    
    
    T_model_arr = data.time_series.T_time_series_4dmat(...
        data.metadata.crater_data.lat_arr == lat,...
        data.metadata.crater_data.long_arr == long,...
        1, :);
    T_model_arr = squeeze(T_model_arr);
    
    clf
    T_model_matrix = data.time_series.T_time_series_4dmat(...
        data.metadata.crater_data.lat_arr == lat,...
        data.metadata.crater_data.long_arr == long,...
        2:end, :);
    T_model_matrix = squeeze(T_model_matrix);
    legend_h_arr = plot(data.time_series.dtm_arr, T_model_matrix, ':g');
    legend_h_arr = legend_h_arr(1);
    hold on
    legend_h_arr(end+1) = plot(data.time_series.dtm_arr, T_model_arr);
%     legend_h_arr(end+1) = scatter(start_dtm_arr, T_diviner_arr, '>');
%     legend_h_arr(end+1) = scatter(end_dtm_arr, T_diviner_arr, '<');
    xlim([data.time_series.dtm_arr(1)-days(30), data.time_series.dtm_arr(end)+days(15)])
    legend_h_arr(end+1) = scatter(convert_jd_to_dtm(jd_arr), T_diviner_arr, 'kx');
%     drawnow
%     legend_h_arr(end+1) = scatter(convert_ltim_to_dtm(ltim_arr, start_dtm_arr, long), T_diviner_arr, 'rx');
    legend(legend_h_arr, {...
        'Model layer temp.',...
        'Model surface temp.',... % 'T_{Diviner} start dtm','T_{Diviner} end dtm',...
        'T_{Diviner} jd'...
%         'T_{Diviner} ltim'...
        })
    xlabel('Datetime')
    ylabel('Temperature (K)')
    title(sprintf('Model & diviner temperatures for single pixel in "%s" at coordinates (%g, %g)', ls_name, lat, long))
    drawnow
end


%% Layer time series
if layer_time_series
    z_arr = data.metadata.z_arr;
    color_arr = NaN(numel(z_arr),3);
    colormap viridis
    colors = colormap;
    clf
    hold on
    for layer_idx = 1:numel(z_arr)
        color_idx = floor(255*z_arr(layer_idx)/z_arr(end))+1;
        plot(data.time_series.dtm_arr, squeeze(data.time_series.T_time_series_4dmat(11,18,layer_idx,:)), 'Color', colors(color_idx, :))
    end
    colorbar
    caxis([z_arr(1), z_arr(end)]);
end

%% Layer depth profile
if layer_depth_profile
    lat_idx = randi(numel(output.crater_data.lat_arr));11;
    long_idx = randi(numel(output.crater_data.long_arr));18;
    z_arr = data.metadata.z_arr*-1;
    clf
    hold on
    plot(squeeze(max(data.time_series.T_time_series_4dmat(lat_idx, long_idx, :, :), [], 4)), z_arr, 'r')
    plot(squeeze(min(data.time_series.T_time_series_4dmat(lat_idx, long_idx, :, :), [], 4)), z_arr, 'b')
    plot(squeeze(data.extremes.Tmax_3dmat(lat_idx, long_idx, :)), z_arr, 'r:')
    plot(squeeze(data.extremes.Tmin_3dmat(lat_idx, long_idx, :)), z_arr, 'b:')
    plot(squeeze(mean(data.time_series.T_time_series_4dmat(lat_idx, long_idx, :, :), 4)), z_arr, 'g')
    title(sprintf('Depth temperature profile for single pixel in "%s" at coordinates (%g, %g)', ls_name, lat_idx, long_idx))
    xlabel('Temperature (K)')
    ylabel('Depth (m)')
    legend('Maximim', 'Minimum', 'Maximim (5y)', 'Minimum (5y)', 'Mean')
    grid on
end

%% Layer depth profile 3d
if profile_3d_old
    lat_idx = randi(24);11;
    long_idx = randi(40);18;
    z_arr = data.metadata.z_arr*-1;
    clf
    hold on
    colors = jet(size(data.time_series.T_time_series_4dmat, 4));
    for time_idx = 1:size(data.time_series.T_time_series_4dmat, 4)
        color_idx = time_idx;
        plot3(squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, :, time_idx)), z_arr, repmat(data.time_series.dtm_arr(time_idx), size(z_arr)), 'Color', colors(color_idx,:))
        if mod(time_idx, 720) == 0
            drawnow
        end
    end
    title(sprintf('(%d, %d)', lat_idx, long_idx))
    xlabel('Layer temperature (K)')
    ylabel('Layer depth (m)')
    zlabel('Time')
end


%% Layer depth profile 3d 2
if profile_3d
    lat_idx = randi(24);11;
    long_idx = randi(40);18;
    z_arr = data.metadata.z_arr*-1;
    clf
    colormap viridis
    [dtm_matrix, z_matrix] = meshgrid(data.time_series.dtm_arr, z_arr);
    surf(z_matrix, dtm_matrix, squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, :, :)))%, z_matrix); %seconds(dtm_matrix-min(dtm_matrix(:)))+1)
    
    shading interp
    title(sprintf('Model layer temperature time series for pixel in "%s" at coordinates: (%d, %d)', ls_name, lat_idx, long_idx))
    xlabel('Depth (m)')
    ylabel('Time')
    zlabel('Temperature (K)')
    step_colormap(spectral(256),0.5,112);
    caxis('auto')
    c_lim = caxis;    
%     diff = min(abs(T_cutoff - c_lim));
%     caxis(T_cutoff + diff*[-1, 1]);
    colorbar
end

%% Depth condition plot
if stable_condition_T
    T_cutoff = 112;
    
    T_3dmat = data.extremes.Tmax_3dmat;
    
    % T_3dmat(T_3dmat > T_cutoff) = NaN;
    % plot_matrix = T_3dmat;
    % elevation_matrix = output.crater_data.elevation_matrix - reshape(data.metadata.z_arr, 1, 1, []);
    % plot_3d_surface(output.crater_data.ew_matrix, output.crater_data.ns_matrix, elevation_matrix, plot_matrix, ls_name);
    z_arr = data.metadata.z_arr;
    plot_matrix = NaN(size(output.crater_data.elevation_matrix));
    alpha_matrix = zeros(size(plot_matrix));
    for lat_idx = 1:numel(output.crater_data.lat_arr)
        for long_idx = 1:numel(output.crater_data.long_arr)
            z_idx_arr = find(T_3dmat(lat_idx, long_idx, :) < T_cutoff);
            if numel(z_idx_arr) > 0
                plot_matrix(lat_idx, long_idx) = z_arr(z_idx_arr(1));
                alpha_matrix(lat_idx, long_idx) = z_arr(z_idx_arr(end))-z_arr(z_idx_arr(1));
                %             if min(z_idx_arr) + numel(z_idx_arr) - numel(z_arr) ~= 1
                %                 alpha_matrix(lat_idx, long_idx) = 0.5;
                %             else
                %                 alpha_matrix(lat_idx, long_idx) = 1;
                %             end
            end
        end
    end
    plot_h = plot_3d_surface(output.crater_data.ew_matrix, output.crater_data.ns_matrix, output.crater_data.elevation_matrix, plot_matrix, ls_name, 'Depth where T_{max} < T_{cutoff} (m)');
    shading faceted
    title(sprintf('Stability depth for "%s" with cutoff temperature = %gK; transparency represents total depth below cutoff temp.', ls_name, T_cutoff))
    alpha(plot_h, alpha_matrix)
end

%% Depth condition plot 2
if stable_condition_z
    T_cutoff = 112;
    T_3dmat = data.extremes.Tmax_3dmat;
    z_arr = data.metadata.z_arr;
    plot_matrix = NaN([size(output.crater_data.elevation_matrix), 3]);
    alpha_matrix = zeros(size(plot_matrix));
    color_matrix = [viridis(6); [0.5 0.5 0.5]; [1 0.9 0.9]];
    % color_matrix = [[1 1 0]; [0.666 1 0]; [0 1 0]; [0 0.666 0.333]; [0 0.333 0.666]; [0 0 1]; [0.333 0 0.666]; [0.666 0 0.333]; [1 0.9 0.9]];
    for lat_idx = 1:numel(output.crater_data.lat_arr)
        for long_idx = 1:numel(output.crater_data.long_arr)
            z_idx_arr = find(T_3dmat(lat_idx, long_idx, :) < T_cutoff);
            if numel(z_idx_arr) > 0
                z_top = z_arr(z_idx_arr(1));
                z_bottom = z_arr(z_idx_arr(end));
                z_diff = z_bottom - z_top;
                if z_top < 0.25
                    if z_bottom > 3
                        if z_top == 0
                            color = color_matrix(1,:);
                        elseif z_top < 0.1
                            color = color_matrix(2,:);
                        else
                            color = color_matrix(3,:);
                        end
                    elseif z_bottom > 1
                        color = color_matrix(4,:);
                    else
                        color = color_matrix(5,:);
                    end
                elseif z_top < 0.5
                    color = color_matrix(6,:);
                else
                    color = color_matrix(7,:);
                end
                
            else
                color = color_matrix(end,:);
            end
            plot_matrix(lat_idx, long_idx, :) = color;
        end
    end
    elevation_matrix = output.crater_data.elevation_matrix-min(output.crater_data.elevation_matrix(:));
    plot_h = plot_3d_surface(output.crater_data.ew_matrix, output.crater_data.ns_matrix, elevation_matrix, plot_matrix);
    shading faceted
    title(sprintf('Stability depth for "%s" with cutoff temperature = %gK', ls_name, T_cutoff))
    colorbar('off')
    hold on
    legend_h_arr = [];
    for idx = 1:size(color_matrix,1)
        legend_h_arr(end+1) = scatter(NaN, NaN, [], 's', 'filled', 'MarkerFaceColor', color_matrix(idx,:));
    end
    scatter3(0,0,interp2(elevation_matrix,size(elevation_matrix,2)/2+0.5,size(elevation_matrix,1)/2+0.5), 'k.')
    legend_h = legend(legend_h_arr, {'surface to >3m','<10cm to >3m', '<25cm to >3m', '<25cm to >1m', '<25cm to <1m', '<50cm to ...', 'Other', 'No stability'});
    % title(legend_h, 'Stable at depths...')
end


%% Plot landing site locations
if plot_ls_locations
    [elevation_matrix, lat_arr, long_arr, r_moon] = read_lola_height_data(16, [-75, -90]);
    elevation_matrix = elevation_matrix - r_moon;
    clf
    axesm('ortho', 'Origin', [-90,0,0], 'MapLatLimit', [-90,-74]);
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 12 6];
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     margin = 0.05;
%     left = outerpos(1) + ti(1) + margin;
%     bottom = outerpos(2) + ti(2) + margin;
%     ax_width = outerpos(3) - ti(1) - ti(3) - 2*margin;
%     ax_height = outerpos(4) - ti(2) - ti(4) - 2*margin;
%     ax.Position = [left bottom ax_width ax_height];
    pcolorm(lat_arr, long_arr, elevation_matrix);
    c = contourcbar;
    c.Label.String = 'Relative elevation (m)';
    title('Landing sites')
    hold on
    gridm on
    mlabel on
    colormap viridis
%     plabel on
    for idx = 1:8
        switch string(idx)
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
        end
        scatterm(ls_lat, ls_long, [], 'xk')
        lat_offset = 0;
        long_offset = 0;
        text_str = sprintf(' %i', idx);
        textm(ls_lat+lat_offset, ls_long+long_offset, 0, text_str, 'Color', [0 0 0])
    end
end


%% Compare deep surface T
if compare_deep_surface_T
    T_deep_arr = double(reshape(mean(data.time_series.T_time_series_4dmat(:, :, end, :), 4), [], 1));
    T_max_arr = double(reshape(max(data.time_series.T_time_series_4dmat(:, :, 1, :), [], 4), [], 1));
    T_min_arr = double(reshape(min(data.time_series.T_time_series_4dmat(:, :, 1, :), [], 4), [], 1));
    T_mean_arr = double(reshape(mean(data.time_series.T_time_series_4dmat(:, :, 1, :), 4), [], 1));


    clf
    hold on
    scatter(T_max_arr, T_deep_arr,  'rx')
    scatter(T_min_arr, T_deep_arr, 'bx')
    scatter(T_mean_arr, T_deep_arr, 'gx')
    lim = [xlim; ylim];
    lim = [max(lim(:,1)), min(lim(:,2))];
    plot(lim, lim, 'k:')
    
    fit_min = fit(T_min_arr, T_deep_arr, 'poly1');
    x_arr = min(lim):0.1:max(lim);
    y_arr = fit_min.p1*x_arr + fit_min.p2;
    
    fprintf('p1 = %g\tp2 = %g\t(%s)\n', fit_min.p1, fit_min.p2, strrep(ls_name, '\_', '_'))
    
%     y_arr(y_arr > max(lim) | y_arr < min(lim)) = NaN;
%     plot(x_arr, y_arr, 'k')
    
    
    xlabel('Surface temperature (K)')
    ylabel('Deep temperature (K)')
    axis equal
    box on
    grid on
    xticks(yticks)
    yticks(xticks)
    legend('Max', 'Min', 'Mean')
    
end