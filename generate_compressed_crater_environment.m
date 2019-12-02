function data = generate_compressed_crater_environment(crater_name, compression_ratio)
% GENERATE_COMPRESSED_CRATER_ENVIRONMENT compresses high latitude crater
% environment so that NS and EW spatial resolution is similar
%
% data = generate_compressed_crater_environment(crater_name) generates data
% file of compressed crater environment for crater_name which is saved and
% returned (uses automatic compression_ratio for equal NS and EW
% resolutions)
%
% data = generate_compressed_crater_environment(crater_name,
% compression_ratio) uses manual compression ratio in EW direction

source_name = create_static_path(sprintf('crater_environments/%s.mat', crater_name));
target_name = create_static_path(sprintf('crater_environments/%s_compressed.mat', crater_name));

data = load(source_name);
data = data.data;
ew_points = size(data.ew_matrix,2);
ns_points = size(data.ns_matrix,1);
ew_length = data.ew_matrix(round(ns_points/2),end)...
    - data.ew_matrix(round(ns_points/2),1);
ns_length = data.ns_matrix(1,1) - data.ns_matrix(end,1);

if nargin < 2
    compression_ratio = round((ew_points/ew_length)/(ns_points/ns_length));
    for i = 0:compression_ratio/3
        % try to ensure compression ratio gives integer divisions
        if mod(ew_points, compression_ratio-i) == 0
            compression_ratio = compression_ratio - i;
            break
        end
    end
end
new_ew_points = ew_points/compression_ratio;
if compression_ratio == 1
    fprintf('No compression necessary')
    return
else
    fprintf('Compressing crater data by %d in ew direction... ', compression_ratio)
end
for lat_idx = 1:numel(data.lat_arr)
    for long_idx = 1:new_ew_points
        old_long_idx_arr = ((long_idx-1)*compression_ratio + 1):long_idx*compression_ratio;
        if lat_idx == 1
            long_arr(long_idx) = mean(data.long_arr(old_long_idx_arr));
        end
        
        ew_matrix(lat_idx, long_idx) = mean(data.ew_matrix(lat_idx, old_long_idx_arr));
        ns_matrix(lat_idx, long_idx) = mean(data.ns_matrix(lat_idx, old_long_idx_arr));
        elevation_matrix(lat_idx, long_idx) = mean(data.elevation_matrix(lat_idx, old_long_idx_arr));
        
        Tmax_matrix(lat_idx, long_idx) = mean(data.Tmax_matrix(lat_idx, old_long_idx_arr));
        Tmin_matrix(lat_idx, long_idx) = mean(data.Tmin_matrix(lat_idx, old_long_idx_arr));
        
        Tmax_ltim_matrix(lat_idx, long_idx) = mode(data.Tmax_ltim_matrix(lat_idx, old_long_idx_arr));
        Tmax_dtm_matrix(lat_idx, long_idx) = mode(data.Tmax_dtm_matrix(lat_idx, old_long_idx_arr));
        Tmin_ltim_matrix(lat_idx, long_idx) = mode(data.Tmin_ltim_matrix(lat_idx, old_long_idx_arr));
        Tmin_dtm_matrix(lat_idx, long_idx) = mode(data.Tmin_dtm_matrix(lat_idx, old_long_idx_arr));
    end
end
data.long_arr = long_arr;
data.ew_matrix = ew_matrix;
data.ns_matrix = ns_matrix;
data.elevation_matrix = elevation_matrix;
data.Tmax_matrix = Tmax_matrix;
data.Tmax_ltim_matrix = Tmax_ltim_matrix;
data.Tmax_dtm_matrix = Tmax_dtm_matrix;
data.Tmin_matrix = Tmin_matrix;
data.Tmin_ltim_matrix = Tmin_ltim_matrix;
data.Tmin_dtm_matrix = Tmin_dtm_matrix;
data.compression_ratio = compression_ratio;
data.compressed_ppd = data.ppd/compression_ratio;
data.compression_created = datetime;
data.compression_description = 'Data for crater, compressed in ew direction to remove unnecessary accuracy at large latitudes';
save(target_name, 'data')
fprintf('Done\n')
end