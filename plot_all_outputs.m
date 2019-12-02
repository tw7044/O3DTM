% PLOT_ALL_OUTPUTS plots all landing site model outputs on same graph

lat_limit = -80;
[plot_matrix, lat_arr, long_arr, r_moon] = read_lola_height_data(16, [lat_limit, -90]);
plot_matrix = plot_matrix - min(plot_matrix(:));
plot_matrix = NaN(size(plot_matrix));
% [plot_matrix, lat_arr, long_arr] = read_diviner_temperatures('max', 16);
% [long_matrix, lat_matrix] = meshgrid(long_arr, lat_arr);
% plot_matrix = plot_matrix(lat_matrix <= lat_limit);
% lat_arr = lat_arr(lat_arr <= lat_limit);
% plot_matrix = reshape(plot_matrix, [numel(lat_arr), numel(long_arr)]);
clf
axesm('ortho', 'Origin', [-90,0,0], 'MapLatLimit', [-90,lat_limit+0.5], 'MLineLocation', 30, 'PLineLocation', 5, 'FontSize', 6, 'LabelRotation', 'on', 'LabelFormat', 'signed');
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
pcolorm(lat_arr, long_arr, plot_matrix);
c = colorbar;
c.Label.String = 'T (K)';
hold on
gridm on
mlabel on
% colormap gray
% caxis([60, 330])
for res = [4 16 128]
    for idx = 1:8
        load(create_static_path(sprintf('outputs/ls_temperatures_extremes/ls%i_%ippd_wide_compressed_temperatures_extremes.mat', idx, res)), 'data');
        
        plot_matrix = data.extremes.Tmax_3dmat(:,:,1);
%         plot_matrix = plot_matrix - data.metadata.crater_data.Tmax_matrix;
        pcolorm(data.metadata.crater_data.lat_arr, data.metadata.crater_data.long_arr, plot_matrix);
        step_colormap(spectral, 0.5, 112);
%         step_colormap(spectral, 0.5, 0);

        data = data.metadata.crater_data;
        [long, lat] = meshgrid(data.long_arr, data.lat_arr);
        lat = reshape(lat, [], 1);
        long = reshape(long, [], 1);
        k = boundary(lat, long);
        plotm(lat(k), long(k), 'w:')
        
        drawnow
    end
end
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
    lat_offset = 0;
    long_offset = 0;
    text_str = sprintf(' %i', idx);
    if idx == 5
        lat_offset = 0.4;
    elseif idx==7
        lat_offset = -0.4;
    end
    scatterm(ls_lat, ls_long, [], 'xk')
    textm(ls_lat+lat_offset, ls_long+long_offset, 0, text_str, 'Color', [0 0 0])
end
drawnow

box off
axis off