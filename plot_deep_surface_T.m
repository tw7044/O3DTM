% PLOT_DEEP_SURFACE_T plots comparison of model deep and surface
% temperatures

clf
hold on
box on
grid on
% axis equal
xlabel('Surface temperature (K)')
ylabel('Minumum stable sub-surface temperature (K)')

% lim = [0, 250];
xlim([0,inf]);
ylim([0,220]);

for idx = 1:8
    for res = [4, 16, 128]
        for ls_name_end = {'', '_wide_compressed'}
            ls_name = sprintf('ls%i_%ippd%s', idx, res, ls_name_end{1});
            plot_subroutine(ls_name)
        end
    end
end

% for idx = 1:8
%     for res = [4, 128]
%         for par = [0.02, 0.1]
%             ls_name = sprintf('ls%i_%ippd_wide_compressed_temperatures_H=%g_extremes', idx, res, par);
%             try
%                 if par == 0.1
%                     plot_subroutine(ls_name, 'x')
%                 else
%                     plot_subroutine(ls_name)
%                 end
%             catch
%                 fprintf('\b XXX \n')
%             end
%         end
%     end
% end

% plot_subroutine('ls1_4ppd_wide_compressed_temperatures_H=0.1_extremes')
% plot_subroutine('ls1_4ppd_wide_compressed_temperatures_H=0.02_extremes', 'x')



function plot_subroutine(ls_name, plot_shape)
if nargin == 1
    plot_shape = '.';
end
fprintf('%s...\n', ls_name)
try
    if contains(ls_name, '=')
        load(create_static_path(sprintf('outputs/ls_temperatures_extremes/parameters/%s.mat', ls_name)), 'data');
    else
        load(create_static_path(sprintf('outputs/ls_temperatures_extremes/%s_temperatures_extremes.mat', ls_name)), 'data');
    end
catch
    fprintf('\b NO FILE\n')
    return
end

if data.metadata.created < datetime(2018, 9, 8)
    fprintf('\b OLD FILE\n')
    return
end

T_stable_arr = reshape(min(data.extremes.Tmax_3dmat(:, :, :), [], 3), [], 1);
T_max_arr = reshape(data.extremes.Tmax_3dmat(:, :, 1), [], 1);
T_min_arr = reshape(data.extremes.Tmin_3dmat(:, :, 1), [], 1);
T_deep_arr = double(reshape(data.time_series_extremes.Tmean_3dmat(:,:,end), [], 1));
T_mean_arr = double(reshape(data.time_series_extremes.Tmean_3dmat(:,:,1), [], 1));

if nargin == 1
    sz = 6;
    scatter(T_max_arr, T_stable_arr, sz, sprintf('%sr', plot_shape))
    scatter(T_min_arr, T_stable_arr, sz, sprintf('%sb', plot_shape))
    scatter(T_mean_arr, T_stable_arr, sz, sprintf('%sg', plot_shape))
%     scatter(T_deep_arr, T_stable_arr, sz, sprintf('%sk', plot_shape))
else
    sz = 6;
    scatter(T_max_arr, T_stable_arr, sz, sprintf('%sk', plot_shape))
    scatter(T_min_arr, T_stable_arr, sz, sprintf('%sc', plot_shape))
    scatter(T_mean_arr, T_stable_arr, sz, sprintf('%sy', plot_shape))
%     scatter(T_deep_arr, T_stable_arr, sz, sprintf('%sb', plot_shape))
end

lim = [xlim; ylim];
lim = [max(lim(:,1)), min(lim(:,2))];
plot(lim, lim, 'Color', 0.8*[1 1 1])

h = refline(0, 112);
h.Color = [0 0 0];
legend('Max surface T', 'Min surface T', 'Mean surface T')%, 'Deep T')
drawnow
end