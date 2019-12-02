% PLOT_PAPER_FIGURES plots figures used in paper

%% CONTROL (set all = 1 to plot all graphs)
enable_plot_multi_profiles = 1;
enable_plot_ls_profile = 0;
enable_plot_ls_accuracy = 0;
enable_plot_ls_stability = 0;
enable_plot_ls_locations = 0;
enable_plot_ls_resolutions = 0;
enable_plot_ls_depth_profile = 0;
enable_plot_3d_profile = 0;
enable_plot_parameter_sensitivity = 0;
enable_plot_min_T_z_histogram = 0;
enable_plot_stable_T_z_histogram = 0;
enable_plot_T_z_heatmap = 0;
enable_plot_stable_T_z_heatmap = 0;
enable_plot_blagg_error = 0;
enable_plot_1d_3d_error = 0;
enable_plot_density_profile = 0;

%% CODE
clf
colormap viridis
prepare_paper_size

%% Plot multi profiles
if enable_plot_multi_profiles
    ww = 8;
    hh = 8;
    prepare_paper_size(hh,ww)
    ppd_list = [128, 16, 4];
    counter = 1;
    ax_list = cell(3,2);
    Tmax = NaN(3,2);
    Tmin = NaN(3,2);
    for ppd_idx = 1:numel(ppd_list)
        ppd = ppd_list(ppd_idx);
        ls_name = sprintf('ls8_%ippd_wide_compressed', ppd);
        data = load_crater_environment(ls_name);
        ew_matrix = data.ew_matrix;
        ns_matrix = data.ns_matrix;
        elevation_matrix = data.elevation_matrix;
        
        T_data = load(create_static_path(sprintf('outputs/ls_temperatures_extremes/%s_temperatures_extremes.mat', ls_name)));
        T_data = T_data.data;
        T_list = {'min', 'max'};
        for T_idx = 1:numel(T_list)
            T_mode = T_list{T_idx};
            bottom_list = [0.6, 0.3, 0.0] + 0.01;
            bottom = bottom_list(ppd_idx);
            left_list = [0.0, 0.5];
            left = left_list(T_idx);
            pos = [left, bottom, 0.5, 0.3];
            ax_list{ppd_idx, T_idx} = subplot('Position', pos);
            
            switch T_mode
                case 'max'
                    plot_matrix = T_data.extremes.Tmax_3dmat(:,:,1);
                case 'min'
                    plot_matrix = T_data.extremes.Tmin_3dmat(:,:,1);
            end
            Tmax(ppd_idx, T_idx) = max(plot_matrix(:));
            Tmin(ppd_idx, T_idx) = min(plot_matrix(:));
            plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, plot_matrix);
            ls_loc_subroutine('r.', ls_name)
            if counter ~= 5
                ylabel('')
            end
            zticks([])
            ticks = [-80, 0, 80];
            if ppd == 128
                ticks = [-2, 0, 2];
            elseif ppd == 16
                ticks = [-40, 0, 40];
            end
            xticks(ticks)
            yticks(ticks)            
            counter = counter + 1;
            
            
            if ppd_idx == 3
                Tmax_ = max(Tmax, [], 1);
                Tmin_ = min(Tmin, [], 1);
                clim = [Tmin_(T_idx), Tmax_(T_idx)];
                for idx = 1:3
                    ax_h = ax_list{idx, T_idx};
                    caxis(ax_h, clim)
                    colormap(ax_h, viridis)
                end
                
                left_list = [0, 0.5] + 0.025;
                left = left_list(T_idx);
                pos = [left, 0.93, 0.45, 0.01];
                c = colorbar('northoutside');
                if T_idx == 1
                    c.Label.String = sprintf('Minimum surface temperature (K)');
                else
                    c.Label.String = sprintf('Maximum surface temperature (K)');
                end
                c.Limits = clim;
                c.Position = pos;
            end
        end
    end

    save_paper_figure('ls_surface_temperature_sumary', hh, ww)
end

%% Plot ls profile
if enable_plot_ls_profile
    ls_name = 'ls8_16ppd_wide_compressed';
    data = load_crater_environment(ls_name);
    ew_matrix = data.ew_matrix;
    ns_matrix = data.ns_matrix;
    elevation_matrix = data.elevation_matrix;
    elevation_matrix = elevation_matrix - min(elevation_matrix(:));
    plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, elevation_matrix/1e3);
    ls_loc_subroutine('r.')
    
    c = colorbar('northoutside');
    c.Label.String = sprintf('Relative elevation (km)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_profile')
end

%% Plot ls accuracy
if enable_plot_ls_accuracy
    clf
    ls_name = 'ls8_16ppd_wide_compressed';
    load(create_static_path(sprintf('outputs/ls_temperatures/diviner_fit/%s_temperatures_diviner_fit.mat', ls_name)), 'data');
    diviner_ls_name = data.metadata.crater_data.crater_name;
    if contains(ls_name, 'compressed')
        diviner_ls_name = sprintf('%s_compressed', diviner_ls_name);
    end
    diviner_data = load(create_static_path(sprintf('crater_environments/diviner_temperatures/%s_diviner_temperatures.mat', diviner_ls_name)));
    diviner_data = diviner_data.output;
    z_arr = data.metadata.z_arr;
    lat_arr = data.metadata.crater_data.lat_arr;
    long_arr = data.metadata.crater_data.long_arr;
    elevation_matrix = data.metadata.crater_data.elevation_matrix;
    elevation_matrix = (elevation_matrix-min(elevation_matrix(:)));
    ew_matrix = data.metadata.crater_data.ew_matrix;
    ns_matrix = data.metadata.crater_data.ns_matrix;
    model_jd_arr = juliandate(data.time_series.dtm_arr);
    
    deal_with_jd_overlap = true;
    
    lat_idx = 33;
    long_idx = 21;
    lat = diviner_data.crater_data.lat_arr(lat_idx);
    long = diviner_data.crater_data.long_arr(long_idx);
    
    T_diviner_arr = diviner_data.T_arr(diviner_data.lat_arr == lat & diviner_data.long_arr == long);
    jd_arr = diviner_data.jd_arr(diviner_data.lat_arr == lat & diviner_data.long_arr == long);
    
    T_model_arr = data.time_series.T_time_series_4dmat(...
        data.metadata.crater_data.lat_arr == lat,...
        data.metadata.crater_data.long_arr == long,...
        1, :);
    T_model_arr = squeeze(T_model_arr);
    
    T_model_matrix = data.time_series.T_time_series_4dmat(...
        data.metadata.crater_data.lat_arr == lat,...
        data.metadata.crater_data.long_arr == long,...
        2:end, :);
    T_model_matrix = squeeze(T_model_matrix);
    legend_h_arr = [];
%     legend_h_arr = plot(data.time_series.dtm_arr, T_model_matrix, 'Color', 0.85*[1 1 1], 'LineStyle', ':');
%     legend_h_arr = legend_h_arr(1);
    hold on
    line_color = lines;
    line_color = line_color(1,:);
    legend_h_arr(end+1) = plot(data.time_series.dtm_arr, T_model_arr, 'Color', line_color);
    xlim([data.time_series.dtm_arr(1), data.time_series.dtm_arr(end)])
    legend_h_arr(end+1) = scatter(convert_jd_to_dtm(jd_arr), T_diviner_arr, 50, 'k.');
    box on
    xlabel('Date')
    ylabel('Surface temperature (K)')
    save_paper_figure('ls_diviner_comparison_time_series')
    
    avg_error_matrix = NaN(size(elevation_matrix));
    abs_error_matrix = NaN(size(elevation_matrix));
    for lat_idx = 1:numel(lat_arr)
        lat = lat_arr(lat_idx);
        for long_idx = 1:numel(long_arr)
            long = long_arr(long_idx);
            diviner_bool_arr = diviner_data.lat_arr == lat & diviner_data.long_arr == long & diviner_data.jd_arr >= min(model_jd_arr) & diviner_data.jd_arr <= max(model_jd_arr);
            diviner_T_arr = diviner_data.T_arr(diviner_bool_arr);
            diviner_jd_arr = diviner_data.jd_arr(diviner_bool_arr);
            avg_err_value = 0;
            abs_err_value = 0;
            T_model_arr = squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, 1, :));
            if deal_with_jd_overlap
                if numel(diviner_T_arr) > 0
                    overlap_start_T_idx = 1;
                    while true
                        jd_value = diviner_jd_arr(overlap_start_T_idx);
                        T_model = interp1(model_jd_arr, T_model_arr, jd_value);
                        overlap_end_T_idx = numel(diviner_T_arr);
                        for T_idx = overlap_start_T_idx:numel(diviner_T_arr)
                            if diviner_jd_arr(T_idx) > jd_value
                                overlap_end_T_idx = T_idx - 1;
                                break
                            end
                        end
                        err_value_arr = T_model - diviner_T_arr(overlap_start_T_idx:overlap_end_T_idx);
                        if max(err_value_arr) > 0 && min(err_value_arr) < 0
                            err_value = 0; % T_model in range of T_diviner
                        else
                            if max(err_value_arr) < 0
                                err_value = max(err_value_arr); % choose err_value closest to 0
                            else
                                err_value = min(err_value_arr);
                            end
                        end
                        avg_err_value = avg_err_value + err_value;
                        abs_err_value = abs_err_value + abs(err_value);
                        
                        overlap_start_T_idx = overlap_end_T_idx + 1;
                        if overlap_start_T_idx > numel(diviner_T_arr)
                            break
                        end
                    end
                end
            else
                for T_idx = 1:numel(diviner_T_arr)
                    T_model = interp1(model_jd_arr, T_model_arr, diviner_jd_arr(T_idx));
                    avg_err_value = avg_err_value + (T_model - diviner_T_arr(T_idx));
                    abs_err_value = abs_err_value + abs(T_model - diviner_T_arr(T_idx));
                end
            end
            if deal_with_jd_overlap
                avg_error_matrix(lat_idx, long_idx) = (avg_err_value/numel(unique(diviner_jd_arr)));
                abs_error_matrix(lat_idx, long_idx) = (abs_err_value/numel(unique(diviner_jd_arr)));
            else
                avg_error_matrix(lat_idx, long_idx) = (avg_err_value/numel(diviner_T_arr));
                abs_error_matrix(lat_idx, long_idx) = (abs_err_value/numel(diviner_T_arr));
            end
        end
    end
    
    clf
    plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, abs_error_matrix);
    label_subroutine
    ls_loc_subroutine
    c = colorbar('northoutside');
    c.Label.String = sprintf('Absolute model temperature error (K)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_diviner_comparison_abs_error_full_colorbar')
    caxis([0,30])
    c.TickLabels(end) = {'>30'};
    save_paper_figure('ls_diviner_comparison_abs_error_tight_colorbar')

    
    clf
    plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, avg_error_matrix);
    label_subroutine
    ls_loc_subroutine
    c = colorbar('northoutside');
    c.Label.String = sprintf('Average model temperature error (K)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    step_colormap(spectral, 1, 0);
    save_paper_figure('ls_diviner_comparison_avg_error_full_colorbar')
    divergent_colormap(30)
    c.TickLabels(end) = {'>30'};
    c.TickLabels(1) = {'<-30'};
    save_paper_figure('ls_diviner_comparison_avg_error_tight_colorbar')

end

%% Plot ls stability
if enable_plot_ls_stability
    ls_name = 'ls8_16ppd_wide_compressed';
    
    colormap viridis
    load(create_static_path(sprintf('crater_environments/ray_tracing/summary/seasons/%s_ray_tracing_summary.mat', ls_name)), 'data');
    ew_matrix = data.ew_matrix;
    ns_matrix = data.ns_matrix;
    elevation_matrix = data.elevation_matrix - min(data.elevation_matrix(:));
    plot_matrix = data.theta_summary.illuminated_time_matrix;
    plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, plot_matrix);
    ls_loc_subroutine('r.')
    c = colorbar('northoutside');
    c.Label.String = sprintf('Fraction of time illuminated');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_stability_map_illuminated_time')
    
    colormap([0.15,0.15,0.15; viridis])
    save_paper_figure('ls_stability_map_illuminated_time_psr')

    
    colormap viridis
    
    ls_name = create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name));
    
    plot_ls_stability(ls_name, 'T');
    colorbar('off')
    label_subroutine
    ls_loc_subroutine('r.')
    c = colorbar('northoutside');
    c.Label.String = sprintf('Minimum sub-surface constant temperature (K)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_stability_map_T_tight_colorbar')
    step_colormap(spectral, 1, 112, true);
    save_paper_figure('ls_stability_map_T_relaxed_colorbar')
    
    plot_ls_stability(ls_name, 'z_min');
    colorbar('off')
    label_subroutine
    ls_loc_subroutine('r.')
    c = colorbar('northoutside');
    c.Label.String = sprintf('Depth of minimum sub-surface constant temperature (m)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_stability_map_z_min_full_color_range')
    caxis('auto')
    save_paper_figure('ls_stability_map_z_min_auto_color_range')
    
    plot_ls_stability(ls_name, 'z_top');
    colorbar('off')
    label_subroutine
    ls_loc_subroutine('r.')
    c = colorbar('northoutside');
    c.Label.String = sprintf('Minimum stable ice depth (m)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_stability_map_z_top_full_color_range')
    caxis('auto')
    save_paper_figure('ls_stability_map_z_top_auto_color_range')

    
    plot_ls_stability(ls_name, 'z_bottom');
    colorbar('off')
    label_subroutine
    ls_loc_subroutine('r.')
    c = colorbar('northoutside');
    c.Label.String = sprintf('Maximum stable ice depth (m)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_stability_map_z_bottom_full_color_range')
    caxis('auto')
    save_paper_figure('ls_stability_map_z_bottom_auto_color_range')

    plot_ls_stability(ls_name, 'z_range');
    colorbar('off')
    label_subroutine
    ls_loc_subroutine('r.')
    c = colorbar('northoutside');
    c.Label.String = sprintf('Total stable ice depth range (m)');
    c.Position = [0.05, 0.84, 0.9, 0.03];
    ax = gca;
    ax_position = ax.Position;
    ax_position(2) = 0.03;
    ax.Position = ax_position;
    save_paper_figure('ls_stability_map_z_range_full_color_range')
    caxis('auto')
    save_paper_figure('ls_stability_map_z_range_auto_color_range')
end


%% Plot landing site locations
if enable_plot_ls_locations
    lat_limit = -75;
    [elevation_matrix, lat_arr, long_arr, r_moon] = read_lola_height_data(16, [lat_limit, -90]);
    elevation_matrix = elevation_matrix - min(elevation_matrix(:));
    clf
    axesm('ortho', 'Origin', [-90,0,0], 'MapLatLimit', [-90,lat_limit+0.5], 'MLineLocation', 45, 'FontSize', 6, 'LabelRotation', 'on', 'LabelFormat', 'signed');
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 12 6];
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    margin = 0.02;
    left = outerpos(1) + ti(1) + margin;
    bottom = outerpos(2) + ti(2) + margin;
    ax_width = outerpos(3) - ti(1) - ti(3) - 2*margin;
    ax_height = outerpos(4) - ti(2) - ti(4) - 2*margin;
    ax.Position = [left bottom ax_width ax_height];
    pcolorm(lat_arr, long_arr, elevation_matrix/1e3);
    c = colorbar;
    c.Label.String = 'Relative elevation (km)';
    hold on
    gridm on
    mlabel on
    colormap viridis
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
        if idx == 5
            lat_offset = 0.4;
        elseif idx==7
            lat_offset = -0.4;
        end
        textm(ls_lat+lat_offset, ls_long+long_offset, 0, text_str, 'Color', [0 0 0])
        
        
%         for res = [4 16 128]
%             data = load_crater_environment(sprintf('ls%i_%ippd_wide_compressed', idx, res));
%             [long, lat] = meshgrid(data.long_arr, data.lat_arr);
%             lat = reshape(lat, [], 1);
%             long = reshape(long, [], 1);
%             k = boundary(lat, long);
%             if res == 4
%                 plotm(lat(k), long(k), 'r')
%             elseif res == 16
%                 plotm(lat(k), long(k), 'w')
%             else
%                 plotm(lat(k), long(k), 'b')
%             end
%         end
    end
    box off
    axis off
    save_paper_figure('landing_site_map')
end


%% PLot landing site resolutions
if enable_plot_ls_resolutions
    ls_name = 'ls8';
    colormap viridis
    paper_height = 2.6;
    prepare_paper_size(paper_height+0.5)
    for stage = [1 2]
        for res = [4, -4, 16, 128, -128]
            if res > 0
                data = load_crater_environment(sprintf('%s_%ippd_wide_compressed', ls_name, res));
            else
                data = load_crater_environment(sprintf('%s_%ippd', ls_name, -res));
            end
            if res == 4 && stage == 1
                clf
                hold on
                height_offset = min(data.elevation_matrix(:));
            end
            x = data.ew_matrix/1e3;
            y = data.ns_matrix/1e3;
            if stage == 1
                pcolor(x, y, (data.elevation_matrix - height_offset)/1e3)
                shading flat
            elseif stage == 2
                x = reshape(x, [], 1);
                y = reshape(y, [], 1);
                k = boundary(x, y);
                if res > 0
                    plot(x(k), y(k), 'r')
                else
                    plot(x(k), y(k), 'w:')
                end
            end
            if res == 4 && stage == 1
                axis equal
                box on
                lim = caxis;
                c = colorbar;
                c.Label.String = 'Relative elevation (km)';
                xlabel('East (km) \rightarrow')
                ylabel('North (km) \rightarrow')
            else
                caxis(lim)
            end
        end
    end
    yticks(xticks)
    save_paper_figure('landing_site_resolutions', paper_height)
end


%% Plot landing site depth profile
if enable_plot_ls_depth_profile
   ls_name = 'ls8_16ppd_wide_compressed';
   load(create_static_path(sprintf('outputs/ls_temperatures_extremes/%s_temperatures_extremes.mat', ls_name)), 'data');
   
   lat_arr = data.metadata.crater_data.lat_arr;
   long_arr = data.metadata.crater_data.long_arr;
   extremes = data.time_series_extremes;
   z_arr = data.metadata.z_arr;
   
   clf
   hold on
   box on
   
   colors = lines(5);
   
   lat_idx = 8;
   long_idx = 22;
   style = '-';
   T_max_arr = squeeze(extremes.Tmax_3dmat(lat_idx, long_idx, :));
   T_min_arr = squeeze(extremes.Tmin_3dmat(lat_idx, long_idx, :));
   T_mean_arr = squeeze(extremes.Tmean_3dmat(lat_idx, long_idx, :));
   plot(T_min_arr, z_arr, 'Color', colors(1, :), 'LineWidth', 1, 'LineStyle', style)
   plot(T_max_arr, z_arr, 'Color', colors(2, :), 'LineWidth', 1, 'LineStyle', style)
   plot(T_mean_arr, z_arr, 'Color', colors(5, :), 'LineWidth', 1, 'LineStyle', style)
   Tmin = min(T_max_arr);
   fprintf('lat = %g\t long = %g\t style = "%s"\t Tmin = %gK\t stable = %gm\n', lat_arr(lat_idx), long_arr(long_idx), style, Tmin, interp1(T_max_arr, z_arr, 112, 'pchip'))
   
   h = area([min(xlim),112], [0 0], max(z_arr), 'LineStyle', 'none', 'FaceColor', 0.5*[1 1 1]);
   alpha(h, 0.15)
   lat_idx = 26;
   long_idx = 20;
   style = '--';
   T_max_arr = squeeze(extremes.Tmax_3dmat(lat_idx, long_idx, :));
   T_min_arr = squeeze(extremes.Tmin_3dmat(lat_idx, long_idx, :));
   T_mean_arr = squeeze(extremes.Tmean_3dmat(lat_idx, long_idx, :));
   plot(T_min_arr, z_arr, 'Color', colors(1, :), 'LineWidth', 1, 'LineStyle', style)
   plot(T_max_arr, z_arr, 'Color', colors(2, :), 'LineWidth', 1, 'LineStyle', style)
   plot(T_mean_arr, z_arr, 'Color', colors(5, :), 'LineWidth', 1, 'LineStyle', style)
   [Tmin, ii] = min(T_max_arr);
   fprintf('lat = %g\t long = %g\t style = "%s\t Tmin = %gK\t stable = ', lat_arr(lat_idx), long_arr(long_idx), style, Tmin)
   fprintf('%gm to ', interp1(T_max_arr(1:ii), z_arr(1:ii), 112, 'pchip'))
   fprintf('%gm\n', interp1(T_max_arr(ii:end), z_arr(ii:end), 112, 'pchip'))
   
   lat_idx = 16;
   long_idx = 16;
   style = ':';
   T_max_arr = squeeze(extremes.Tmax_3dmat(lat_idx, long_idx, :));
   T_min_arr = squeeze(extremes.Tmin_3dmat(lat_idx, long_idx, :));
   T_mean_arr = squeeze(extremes.Tmean_3dmat(lat_idx, long_idx, :));
   plot(T_min_arr, z_arr, 'Color', colors(1, :), 'LineWidth', 1, 'LineStyle', style)
   plot(T_max_arr, z_arr, 'Color', colors(2, :), 'LineWidth', 1, 'LineStyle', style)
   plot(T_mean_arr, z_arr, 'Color', colors(5, :), 'LineWidth', 1, 'LineStyle', style)
   Tmin = min(T_max_arr);
   fprintf('lat = %g\t long = %g\t style = "%s"\t Tmin = %gK\n', lat_arr(lat_idx), long_arr(long_idx), style, Tmin)
   

   
   xlim([-inf, max(T_max_arr)])
   
   leg = legend('Minimum', 'Maximum', 'Mean');
   leg.Location = 'southeast';
   title(leg, 'Temperature')
   
   xlabel('Temperature (K)')
   ylabel('Depth (m)')
   
   ytickformat('%0.1f')
   ylim([min(z_arr), max(z_arr)])
   ax = gca;
   ax.YDir = 'reverse';
   
   save_paper_figure('T_z_graph')
end


%% Plot 3d temperature profile
if enable_plot_3d_profile
    ls_name = 'ls1_4ppd_wide_compressed';
    load(create_static_path(sprintf('outputs/ls_temperatures/%s_temperatures.mat', ls_name)), 'data');
    
    lat_idx = 16;
    long_idx = 16;
    z_arr = data.metadata.z_arr*-1;
    
    clf
    
    prepare_paper_size(6,8)
    
    colormap viridis
    [dtm_matrix, z_matrix] = meshgrid(data.time_series.dtm_arr, z_arr);
    surf(z_matrix, dtm_matrix, squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, :, :)))%, z_matrix); %seconds(dtm_matrix-min(dtm_matrix(:)))+1)
    
    xlim([min(z_arr), max(z_arr)])
    
    shading interp
    xlabel('Depth (m)')
    ylabel('Time')
    zlabel('Temperature (K)')
    
    view(75,15)
    
    save_paper_figure('3d_temperature_profile', 6, 8)
end



%% Plot parameter sensitivity
if enable_plot_parameter_sensitivity
    global t_wait n_pts lat max_range
    lat = 80;
    t_wait = 'seasons';
    n_pts = 100;
    max_range = 0;
    H_arr = linspace(1e-9, 0.30, n_pts);
    A0_arr = linspace(0, 0.4, n_pts);
    height = 2.1;
    
    parameter_sensitivity_subroutine('A0', A0_arr, 'Bond albedo, A_0')
    save_paper_figure('A0_sensitivity', height, 4)
    
    parameter_sensitivity_subroutine('H', H_arr, 'H parameter (cm)')
    save_paper_figure('H_sensitivity', height, 4)
    
%     parameter_sensitivity_subroutine('A0', A0_arr, 'Bond albedo, A_0')
%     save_paper_figure('A0_sensitivity', 2, 4)
end



%% Plot min T z histogram
if enable_plot_min_T_z_histogram
    ls_name = 'ls8_16ppd_wide_compressed';
    load(create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name)), 'data');

    T_3dmat = data.extremes.Tmax_3dmat;
    z_arr = data.metadata.z_arr;
    lat_arr = data.metadata.crater_data.lat_arr;
    long_arr = data.metadata.crater_data.long_arr;
    elevation_matrix = data.metadata.crater_data.elevation_matrix;
    elevation_matrix = (elevation_matrix-min(elevation_matrix(:)))/1e3;
    ew_matrix = data.metadata.crater_data.ew_matrix/1e3;
    ns_matrix = data.metadata.crater_data.ns_matrix/1e3;
    plot_matrix = NaN(size(elevation_matrix));
    T_matrix = NaN(size(elevation_matrix));
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            dz = 0.01; % 1cm
            z_interp_arr = min(z_arr):dz:max(z_arr);
            T_arr = squeeze(T_3dmat(lat_idx, long_idx, :));
            T_interp_arr = interp1(z_arr, T_arr, z_interp_arr, 'spline');
            [T, value_idx] = min(T_interp_arr);
            T_matrix(lat_idx, long_idx) = T;
            plot_matrix(lat_idx, long_idx) = z_interp_arr(value_idx);
        end
    end
    bin_width = 0.05;
    bin_edges = 0:bin_width:(max(plot_matrix(:))+bin_width);
    clf
    histogram(plot_matrix, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
%     histogram(plot_matrix(T_matrix < 112), bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
%     hold on
%     histogram(plot_matrix(T_matrix >= 112), bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
    xlim([0,inf])
    xlabel('Depth of minimum constant temperature (m)')
    ylabel('Number of surface elements')
%     leg = legend('T < 112K', 'T > 112K');
%     leg.Location = 'northwest';
%     title(leg, 'Min. constant temp.')
    save_paper_figure('z_Tmin_histogram')
end


%% Plot stable T z histogram
if enable_plot_stable_T_z_histogram
    ls_name = 'ls8_16ppd_wide_compressed';
    load(create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name)), 'data');
    
    rt_data = load(create_static_path(sprintf('crater_environments/ray_tracing/summary/seasons/%s_ray_tracing_summary.mat', ls_name)));
    illuminated_time_matrix = rt_data.data.theta_summary.illuminated_time_matrix;
    
    T_3dmat = data.extremes.Tmax_3dmat;
    z_arr = data.metadata.z_arr;
    lat_arr = data.metadata.crater_data.lat_arr;
    long_arr = data.metadata.crater_data.long_arr;
    elevation_matrix = data.metadata.crater_data.elevation_matrix;
    elevation_matrix = (elevation_matrix-min(elevation_matrix(:)))/1e3;
    ew_matrix = data.metadata.crater_data.ew_matrix/1e3;
    ns_matrix = data.metadata.crater_data.ns_matrix/1e3;
    plot_matrix = NaN(size(elevation_matrix));
    T_matrix = NaN(size(elevation_matrix));
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            dz = 0.01; % 1cm
            z_interp_arr = min(z_arr):dz:max(z_arr);
            T_arr = squeeze(T_3dmat(lat_idx, long_idx, :));
            T_interp_arr = interp1(z_arr, T_arr, z_interp_arr, 'spline');
            z_stable_arr = z_interp_arr(T_interp_arr < 112);
            if numel(z_stable_arr) > 0
                plot_matrix(lat_idx, long_idx) = z_stable_arr(1);
            end
        end
    end
    bin_width = 0.02;
    bin_edges = 0:bin_width:(max(plot_matrix(:))+bin_width);
    clf
    histogram(plot_matrix(illuminated_time_matrix == 0), bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
    hold on
    histogram(plot_matrix(illuminated_time_matrix > 0), bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
    xlim([0,inf])
    xlabel('Minimum stable ice depth (m)')
    ylabel('Number of surface elements')
    leg = legend('Permanently shaded', 'Illuminated');
%     leg.Location = 'northwest';
%     title(leg, 'Min. constant temp.')
    save_paper_figure('stable_z_top_histogram')
    
    clf
    bin_width = 0.025;
    bin_edges = 0:bin_width:(max(illuminated_time_matrix(:))+bin_width);
    histogram(illuminated_time_matrix(~isnan(plot_matrix)), bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
    hold on
    histogram(illuminated_time_matrix(isnan(plot_matrix)), bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 1)
    xlim([0,inf])
    xlabel('Fraction of time illuminated')
    ylabel('Number of surface elements')
    leg = legend('T < 112K', 'T > 112K');
    leg.Location = 'northwest';
    title(leg, 'Min. constant temp.')
    save_paper_figure('stable_T_illuminated_histogram')
end


%% Plot stable T z heatmap
if enable_plot_T_z_heatmap
    ls_name = 'ls8_16ppd_wide_compressed';
    load(create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name)), 'data');

    T_3dmat = data.extremes.Tmax_3dmat;
    z_arr = data.metadata.z_arr;
    lat_arr = data.metadata.crater_data.lat_arr;
    long_arr = data.metadata.crater_data.long_arr;
    elevation_matrix = data.metadata.crater_data.elevation_matrix;
    elevation_matrix = (elevation_matrix-min(elevation_matrix(:)))/1e3;
    ew_matrix = data.metadata.crater_data.ew_matrix/1e3;
    ns_matrix = data.metadata.crater_data.ns_matrix/1e3;
    plot_matrix = NaN(size(elevation_matrix));
    T_matrix = NaN(size(elevation_matrix));
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            dz = 0.01; % 1cm
            z_interp_arr = min(z_arr):dz:max(z_arr);
            T_arr = squeeze(T_3dmat(lat_idx, long_idx, :));
            T_interp_arr = interp1(z_arr, T_arr, z_interp_arr, 'spline');
            [T, value_idx] = min(T_interp_arr);
            T_matrix(lat_idx, long_idx) = T;
            plot_matrix(lat_idx, long_idx) = z_interp_arr(value_idx);
        end
    end
    
    colormap(log_colormap(flipud(viridis), 1))
    xbins = min(T_matrix):5:max(T_matrix);
    ybins = 0:0.05:max(plot_matrix(:));
    
    histogram2(T_matrix, plot_matrix, xbins, ybins, 'DisplayStyle', 'tile', 'EdgeColor', 'none');
    h1 = gca;
    set(h1, 'Ydir', 'reverse')
    xlabel('Minimum constant temperature (K)')
    ylabel('Depth of min. const. temperature (m)')
    c = colorbar;
    c.Label.String = 'Number of surface elements';
    grid off
    c = caxis;
    c(1) = 1;
    caxis(c)
    save_paper_figure('z_Tmin_2dhistogram')
end


%% Plot stable T z heatmap
if enable_plot_stable_T_z_heatmap
    ls_name = 'ls8_16ppd_wide_compressed';
    load(create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name)), 'data');
    
    rt_data = load(create_static_path(sprintf('crater_environments/ray_tracing/summary/seasons/%s_ray_tracing_summary.mat', ls_name)));
    illuminated_time_matrix = rt_data.data.theta_summary.illuminated_time_matrix;
    
    T_3dmat = data.extremes.Tmax_3dmat;
    z_arr = data.metadata.z_arr;
    lat_arr = data.metadata.crater_data.lat_arr;
    long_arr = data.metadata.crater_data.long_arr;
    elevation_matrix = data.metadata.crater_data.elevation_matrix;
    elevation_matrix = (elevation_matrix-min(elevation_matrix(:)))/1e3;
    ew_matrix = data.metadata.crater_data.ew_matrix/1e3;
    ns_matrix = data.metadata.crater_data.ns_matrix/1e3;
    plot_matrix = NaN(size(elevation_matrix));
    T_matrix = NaN(size(elevation_matrix));
    for lat_idx = 1:numel(lat_arr)
        for long_idx = 1:numel(long_arr)
            dz = 0.01; % 1cm
            z_interp_arr = min(z_arr):dz:max(z_arr);
            T_arr = squeeze(T_3dmat(lat_idx, long_idx, :));
            T_interp_arr = interp1(z_arr, T_arr, z_interp_arr, 'spline');
            z_stable_arr = z_interp_arr(T_interp_arr < 112);
            if numel(z_stable_arr) > 0
                plot_matrix(lat_idx, long_idx) = z_stable_arr(1);
            end
        end
    end
    
    colormap(log_colormap(flipud(viridis), 1))
    xbins = 0:0.05:max(illuminated_time_matrix(~isnan(plot_matrix)));
    
    dy = 0.025;
    yy = unique(plot_matrix(~isnan(plot_matrix)));
    y0 = -dy*0.75;
    ybins = y0:dy:max(plot_matrix(:));

    histogram2(illuminated_time_matrix, plot_matrix, xbins, ybins, 'DisplayStyle', 'tile', 'EdgeColor', 'none');
    h1 = gca;
    set(h1, 'Ydir', 'reverse')
    xlabel('Fraction of time illuminated')
    ylabel('Minimum stable ice depth (m)')
    c = colorbar;
    c.Label.String = 'Number of surface elements';
    grid off
    c = caxis;
    c(1) = 1;
    caxis(c)
    save_paper_figure('stable_T_illuminated_2dhistogram')
end




%% Plot blagg error
if enable_plot_blagg_error
    crater_data = load('outputs/3d/1d_3d_comparison/blagg_match_max.mat');
    crater_data = crater_data.data;
    
    c_min = min(min(min(crater_data.hayne_1d_3dmat(:,:,1) - crater_data.Tmax_matrix)),...
        min(min(crater_data.full_3d_3dmat(:,:,1) - crater_data.Tmax_matrix)));
    c_max = max(max(max(crater_data.hayne_1d_3dmat(:,:,1) - crater_data.Tmax_matrix)),...
        max(max(crater_data.full_3d_3dmat(:,:,1) - crater_data.Tmax_matrix)));
    c_min = -15;
    c_lim = [c_min, c_max];
    
    local_colors = viridis(256);
    color_extend = round(256*c_max/(-c_min));
    start_color = 1-0.5*(1-local_colors(end, :));
    end_color = [1 1 1];
    start_color = local_colors(end, :);
    for idx = 1:color_extend
        frac = idx/color_extend;
        color = start_color - frac*(start_color - end_color);
        local_colors(end+1, :) = color;
    end
    
    local_colors = spectral(256);
    cutoff = 128 + round(128*c_max/(-c_min));
    local_colors = local_colors(1:cutoff, :);
    colormap(local_colors)
    
    plot_matrix = crater_data.hayne_1d_3dmat(:,:,1) - crater_data.Tmax_matrix;
    real_crater_subroutine(crater_data, plot_matrix)
    text(-2, 8, '\bfA  \rm1D model error')
    zticks([0,1])
    caxis(c_lim)
    colormap(local_colors)
    save_paper_figure('blagg_error_1d', 2.5)
    
    plot_matrix = crater_data.full_3d_3dmat(:,:,1) - crater_data.Tmax_matrix;
    real_crater_subroutine(crater_data, plot_matrix)
    text(-2, 8, '\bfB  \rm3D model error')
    zticks([0,1])
    caxis(c_lim)
    colormap(local_colors)
    save_paper_figure('blagg_error_3d', 2.5)
    
    clf
    axis off
    c = colorbar('north');
    caxis(c_lim)
    c.Label.String = 'Model maximum temperature error, T_{Model} - T_{Diviner} (K)';
    c.Position = [0.05, 0.05, 0.9, 0.15];
    save_paper_figure('blagg_error_colorbar', 0.6, 4)
    
    
    clf
    hold on
    T_1d_matrix = crater_data.hayne_1d_3dmat(:,:,1);
    T_3d_matrix = crater_data.full_3d_3dmat(:,:,1);
    color_arr = flipud(viridis(256));
    diff_arr = linspace(0,max(max(T_3d_matrix - T_1d_matrix)), 256);
    Tdiviner = crater_data.Tmax_matrix;
    for idx = 1:numel(crater_data.Tmax_matrix)
        Td = Tdiviner(idx);
        T1d = T_1d_matrix(idx);
        T3d = T_3d_matrix(idx);
        color = interp1(diff_arr, color_arr, T3d - T1d);
        if isnan(color)
            color = color_arr(1,:);
        end
        plot([Td, Td], [T1d, T3d]-Td, 'Color', color)
        scatter(Td, T3d-Td, 10, color, '^', 'filled')
    end
    text(374.5, 6, '\bfC')
    xlim([0.999*min(Tdiviner(:)), 1.001*max(Tdiviner(:))]);
    plot(xlim, [0,0], 'k--')
    xlabel('Diviner temperature (K)')
    ylabel('T_{Model} - T_{Diviner} (K)')
    box on
    save_paper_figure('blagg_error_graph', 2)
end



%% Plot ls 1d 3d error
if enable_plot_1d_3d_error
    data_3d = load(create_static_path('outputs/ls_temperatures_extremes/diviner_fit/ls8_16ppd_wide_compressed_temperatures_diviner_fit_extremes.mat'));
    data_3d = data_3d.data;
    
    data_1d = load(create_static_path('outputs/ls_1d_temperatures/ls8_16ppd_wide_compressed_1d_temperatures.mat'));
    data_1d = data_1d.data;
    
    T_3d_matrix = squeeze(data_3d.extremes.Tmax_3dmat(:,:,1));
    T_1d_matrix = data_1d.Tmax_matrix;
    T_diviner_matrix = crater_data.Tmax_matrix;
    
    crater_data = data_1d.crater_data;
    
    c_min = min(min(min(T_1d_matrix - T_diviner_matrix)),...
        min(min(T_3d_matrix - T_diviner_matrix)));
    c_max = max(max(max(T_1d_matrix - T_diviner_matrix)),...
        max(max(T_3d_matrix - T_diviner_matrix)));
    c_min = -50;
    c_lim = [c_min, c_max];
    
    caxis(c_lim);
    local_colors = step_colormap(spectral, 1, 0); 
%     local_colors = spectral(256);
%     cutoff = 128 + round(128*c_max/(-c_min));
%     local_colors = local_colors(1:cutoff, :);
    colormap(local_colors)
    
    text_x = -25;
    text_y = 125;
    
    plot_matrix = T_1d_matrix - T_diviner_matrix;
    plot_subroutine(crater_data.ew_matrix, crater_data.ns_matrix, crater_data.elevation_matrix, plot_matrix);
    text(text_x, text_y, '\bfA  \rm1D model error')
    caxis(c_lim)
    colormap(local_colors)
    save_paper_figure('ls_error_1d', 2.5)
    
    plot_matrix = T_3d_matrix - T_diviner_matrix;
    plot_subroutine(crater_data.ew_matrix, crater_data.ns_matrix, crater_data.elevation_matrix, plot_matrix);
    text(text_x, text_y, '\bfB  \rm3D model error')
    caxis(c_lim)
    colormap(local_colors)
    save_paper_figure('ls_error_3d', 2.5)
    
    clf
    axis off
    c = colorbar('north');
    caxis(c_lim)
    c.Label.String = 'Model maximum temperature error, T_{Model} - T_{Diviner} (K)';
    c.Position = [0.05, 0.05, 0.9, 0.15];
    save_paper_figure('ls_error_colorbar', 0.6, 4)
    
    
    clf
    hold on
    diff_arr = linspace(min(min(T_3d_matrix - T_1d_matrix)),max(max(T_3d_matrix - T_1d_matrix)), 256);
    color_arr = lines(5);
    for idx = 1:numel(T_diviner_matrix)
        Td = T_diviner_matrix(idx);
        T1d = T_1d_matrix(idx);
        T3d = T_3d_matrix(idx);
%         color = interp1(diff_arr, color_arr, T3d - T1d);
        
%         if isnan(color)
%             color = color_arr(1,:);
%         end
%         color = 'k';
        
        color = color_arr(2, :);
        shape_str = '^';
        if T3d - T1d < 0
            shape_str = 'v';
            color = color_arr(1, :);
        end
        
%         if abs(T1d - Td) > abs(T3d - Td)
%             color = color_arr(5,:);
%         else
%             color = color_arr(2,:);
%         end
        
        plot([Td, Td], [T1d, T3d]-Td, 'Color', color)

        scatter(Td, T3d-Td, 10, color, shape_str, 'filled')
%         if mod(idx, 100) == 0
%             drawnow
%         end
    end
    text(100, 175, '\bfC')
    xlim([0.999*min(T_diviner_matrix(:)), 1.001*max(T_diviner_matrix(:))]);
    plot(xlim, [0,0], 'k--')
    xlabel('Diviner temperature (K)')
    ylabel('T_{Model} - T_{Diviner} (K)')
    box on
    save_paper_figure('ls_error_graph', 2)
end


%% Density profile
if enable_plot_density_profile
    H_pts = 6;
    H_min = 0.05;
    H_max = 0.30;
    z_min = 0;
    z_max = 0.5;
    z_pts = 100;
    clf
    par = define_parameters;
    rho_d = par.rho_d;
    rho_s = par.rho_s;
    c_arr = viridis(H_pts);
    hold on
    H_arr = linspace(H_min, H_max, H_pts);
    z_arr = linspace(z_min, z_max, z_pts);
    legend_arr = {};
    for H_idx = 1:H_pts
        H = H_arr(H_idx);
        rho_arr = [];
        for z_idx = 1:z_pts
            z = z_arr(z_idx);
            rho_arr(z_idx) = rho_d - (rho_d - rho_s)*exp(-z/H);
        end
        plot(rho_arr, z_arr*100, 'Color', c_arr(H_idx,:), 'LineWidth', 1)
        legend_arr{end+1} = sprintf('%.f cm',H*100);
    end
    % c = colorbar;
    % caxis([H_min, H_max])
    % c.Label.String = 'H parameter (m)';
    xlabel('Density (kg m^{-3})')
    ylabel('Depth (cm)')
    set(gca,'Ydir','reverse')
    l = legend(legend_arr, 'Location', 'southwest');
    title(l, 'H parameter')
    xlim([rho_s, rho_d])
    box on
    
    save_paper_figure('density_profile')
end

%% SUBROUTINES
function label_subroutine
% xh = xlabel('East \rightarrow');
xticks('auto')
yticks('auto')
yh = ylabel('\leftarrow North');
zlabel('')
xlabel('')
title('')
% set(xh, 'Rotation', 35);
set(yh, 'Rotation', -35);
ztickformat('%gkm')
xtickformat('%gkm')
ytickformat('%gkm')
yticks(xticks)
if numel(zticks) == 0
    zticks([0 10])
end
zticks([min(zticks), max(zticks)])
view([-45,45])
box on
end


function prepare_paper_size(height, width)
if nargin < 1 || isempty(height)
    height = 3;
end
if nargin < 2 || isempty(width)
    width = 4;
end
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 width height];
fig.Units = 'inches';
fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
drawnow
end

function parameter_sensitivity_subroutine(variable_name, variable_arr, x_label)
global t_wait n_pts lat max_range
Tmax_arr = NaN(size(variable_arr));
Tmin_arr = NaN(size(variable_arr));
for idx = 1:round(n_pts/2)
    fprintf('_')
end
fprintf('\n')
for idx = 1:n_pts
    if mod(idx, 2) == 0
        fprintf('\b:')
    else
        fprintf('.')
    end
    Tmax_out = hayne_1d_model(lat, t_wait, 'max', 0, 0, struct(variable_name, variable_arr(idx)));
    Tmax_arr(idx) = Tmax_out(1);
    Tmin_out = hayne_1d_model(lat, t_wait, 'min', 0, 0, struct(variable_name, variable_arr(idx)));
    Tmin_arr(idx) = Tmin_out(1);
end
fprintf('\n')
if strcmp(variable_name, 'H')
    variable_arr = variable_arr*100;
end
clf
hold on
yyaxis left
plot(variable_arr, Tmin_arr, 'LineWidth', 1)
xlabel(x_label)
ylabel('Min. temperature (K)')
if range(Tmin_arr) > max_range
    max_range = ceil(range(Tmin_arr)+1);
end
yyaxis right
plot(variable_arr, Tmax_arr, '--', 'LineWidth', 1)
ylabel('Max. temperature (K)')
if range(Tmax_arr) > max_range
    max_range = ceil(range(Tmax_arr)+1);
end
legend('T_{min}', 'T_{max}')
box on
range_offset = 0.5*[-max_range, max_range];
yyaxis left
ylim(round(min(Tmin_arr) + 0.5*range(Tmin_arr)) + range_offset)
yyaxis right
ylim(round(min(Tmax_arr) + 0.5*range(Tmax_arr)) + range_offset)
end

function real_crater_subroutine(crater_data, plot_matrix)
clf
elevation_matrix = crater_data.elevation_matrix - referenceSphere('moon').Radius;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));
plot_3d_surface(crater_data.ew_matrix/1e3, crater_data.ns_matrix/1e3, elevation_matrix/1e3, squeeze(plot_matrix))
colorbar('off')
grid off
label_subroutine
shading faceted
box on
ax = gca;
ax_position = ax.Position;
ax_position(2) = 0.05;
ax_position(4) = 0.9;
ax.Position = ax_position;
end


function ls_loc_subroutine(style, ls_name)
if nargin == 0
    style = 'r.';
end
if nargin < 2
    ls_name = 'ls8_16ppd_wide_compressed';
end
data = load_crater_environment(ls_name);
ls_lat = -82.7;
ls_long = 33.50;
lat_idx = interp1(data.lat_arr, 1:numel(data.lat_arr), ls_lat);
long_idx = interp1(data.long_arr, 1:numel(data.long_arr), ls_long);
% fprintf('LAT  =%g\t min=%g\t max=%g\t len=%g\t idx=%g\n', ls_lat, min(lat_arr), max(lat_arr), numel(lat_arr), lat_idx)
% fprintf('LONG =%g\t min=%g\t max=%g\t len=%g\t idx=%g\n', ls_long, min(long_arr), max(long_arr), numel(long_arr), long_idx)

elevation_matrix = data.elevation_matrix - min(data.elevation_matrix(:));
elev = interp2(elevation_matrix, long_idx, lat_idx);
ew = interp2(data.ew_matrix, long_idx, lat_idx);
ns = interp2(data.ns_matrix, long_idx, lat_idx);

elev = elev/1e3;
ew = ew/1e3;
ns = ns/1e3;
hold on
scatter3(ew, ns, elev, 25, style)
hold off
end

function plot_h = plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, plot_matrix)
elevation_matrix = elevation_matrix - min(elevation_matrix(:));
plot_h = surf(ew_matrix/1e3, ns_matrix/1e3, elevation_matrix/1e3, plot_matrix, 'FaceColor', 'interp');
plot_h.LineWidth = 0.1;
plot_h.EdgeAlpha = 0.33;
grid off
axis equal
shading faceted
label_subroutine
end
