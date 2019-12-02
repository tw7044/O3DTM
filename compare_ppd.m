% COMPARE_PPD creates a plot comparing diviner temperature values at
% different resolutions

mode_1 = true;
mode_2 = false;



if mode_1
    lat = -82.625;
    long = 33.375;
    buffer = 0.125;
    
    data_4 = load(create_static_path('crater_environments/diviner_temperatures/ls8_4ppd_wide_compressed_diviner_temperatures.mat'));
    data_4 = data_4.output;
    data_128 = load(create_static_path('crater_environments/diviner_temperatures/ls8_128ppd_wide_compressed_diviner_temperatures.mat'));
    data_128 = data_128.output;
    clf
    hold on
    
    long_4 = data_4.uncompressed_long_arr;
    lat_4 = data_4.lat_arr;
    dtm_4 = convert_jd_to_dtm(data_4.jd_arr);
    T_4 = data_4.T_arr;
    
    long_128 = data_128.uncompressed_long_arr;
    lat_128 = data_128.lat_arr;
    dtm_128 = convert_jd_to_dtm(data_128.jd_arr);
    T_128 = data_128.T_arr;

    bool_4 = abs(long_4 - long) < buffer & abs(lat_4 - lat) < buffer & dtm_4 < max(dtm_128);
    long_4 = long_4(bool_4);
    lat_4 = lat_4(bool_4);
    dtm_4 = dtm_4(bool_4);
    T_4 = T_4(bool_4);
    
    bool_128 = abs(long_128 - long) < buffer & abs(lat_128 - lat) < buffer;
    long_128 = long_128(bool_128);
    lat_128 = lat_128(bool_128);
    dtm_128 = dtm_128(bool_128);
    T_128 = T_128(bool_128);
    
    scatter(dtm_128, T_128, 1, '.k')
   
    T_mean = NaN(size(T_4));
    T_other = NaN(size(T_4));
    for idx = 1:numel(T_4)
        dtm = dtm_4(idx);
        bool_arr = abs(dtm_128 - dtm) < days(1);
        T_mean(idx) = mean(T_128(bool_arr));
        T_other(idx) = (sum(T_128(bool_arr).^4)/sum(bool_arr)).^0.25;
    end
    scatter(dtm_4, T_mean, 'ok')
    
    scatter(dtm_4, T_4, 'xr')
    scatter(dtm_4, T_other, 'og')

    legend('128ppd', '128ppd mean', '4ppd', 'avg flux')
    xlabel('Date')
    ylabel('Temperature')
    zlabel('Latitude')
    title(sprintf('Diviner measurements for (%g, %g)', lat, long))
end


if mode_2
    lat_min = 0.125;
    lat_max = 0.625;
    long_min = 0.125;
    long_max = 0.625;
    
    % fprintf('Loading elevation data [4 ...')
    % % [elevation_matrix_4, lat_arr_4, long_arr_4] = read_lola_height_data(4, [lat_min, lat_max], [long_min, long_max]);
    % fprintf('] [128 ...')
    % % [elevation_matrix_128, lat_arr_128, long_arr_128] = read_lola_height_data(128, [lat_min, lat_max], [long_min, long_max]);
    % fprintf('] Done\n')
    % %
    % % [Tmax_matrix_4, Tmin_matrix_4, Tmax_ltim_matrix_4, Tmin_ltim_matrix_4, Tmax_dtm_matrix_4, Tmin_dtm_matrix_4] = find_diviner_extreme_temperatures(4, lat_arr_4, long_arr_4);
    % % [Tmax_matrix_128, Tmin_matrix_128, Tmax_ltim_matrix_128, Tmin_ltim_matrix_128, Tmax_dtm_matrix_128, Tmin_dtm_matrix_128] = find_diviner_extreme_temperatures(128, lat_arr_128, long_arr_128);
    load('compare_ppd_environment')
    
    clf
    colormap viridis
    
    [long_matrix_4, lat_matrix_4] = meshgrid(long_arr_4, lat_arr_4);
    % subplot(1,3,1)
    % pcolor(long_matrix_4, lat_matrix_4, Tmax_matrix_4);
    % ylim([lat_min, lat_max])
    % xlim([long_min, long_max])
    % shading flat
    % hold on
    % title('4ppd')
    value_matrix_4 = Tmax_matrix_4;
    value_matrix_128 = Tmax_matrix_128;
    
    color_limits = [min([value_matrix_4(:); value_matrix_128(:)]), max([value_matrix_4(:); value_matrix_128(:)])];
    
    subplot(1,3,1)
    [long_matrix_128, lat_matrix_128] = meshgrid(long_arr_128, lat_arr_128);
    value_matrix_interp = interp2(long_matrix_4, lat_matrix_4, value_matrix_4, long_matrix_128, lat_matrix_128);
    pcolor(long_matrix_128, lat_matrix_128, value_matrix_interp);
    shading flat
    ylim([lat_min, lat_max])
    xlim([long_min, long_max])
    hold on
    scatter(reshape(long_matrix_4,1,[]),reshape(lat_matrix_4,1,[]),[],reshape(value_matrix_4,1,[]), 'filled', 'MarkerEdgeColor', 'k')
    title('4ppd data (circles) & interpolated values to 128ppd')
    caxis(color_limits)
    colorbar('location', 'southoutside')
    
    subplot(1,3,2)
    pcolor(long_matrix_128, lat_matrix_128, value_matrix_128);
    shading flat
    ylim([lat_min, lat_max])
    xlim([long_min, long_max])
    title('128ppd data')
    caxis(color_limits)
    colorbar('location', 'southoutside')
    
    ax3 = subplot(1,3,3);
    pcolor(long_matrix_128, lat_matrix_128, value_matrix_interp - value_matrix_128);
    shading flat
    ylim([lat_min, lat_max])
    xlim([long_min, long_max])
    title('Interpolated - 128ppd data')
    colorbar('location', 'southoutside')
    caxis('auto')
    lim = max(abs(caxis));
    caxis(lim*[-1 1])
    colormap(ax3, spectral)
end