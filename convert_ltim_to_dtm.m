function output_dtm = convert_ltim_to_dtm(ltim_input, start_dtm_input, long_input)
% CONVERT_LTIM_TO_DTM converts a local lunar time & diviner start time to
% datetime (UTC).
%
% output_dtm = CONVERT_LTIM_TO_DTM(ltim_input, start_dtm_input, long_input)
% converts local time matrix, start datetime matrix and long_matrix to
% datetime of measurements.

[horizons_sun_dtm_arr, horizons_sun_long_arr] = read_horizons_data();
output_dtm = NaT(size(ltim_input));

if numel(long_input) == 1 && numel(ltim_input) > 1
    long_input = long_input*ones(size(ltim_input));
end
for ltim_idx = 1:numel(ltim_input)
    ltim = ltim_input(ltim_idx);
    start_dtm = start_dtm_input(ltim_idx);
    long = long_input(ltim_idx);
    
    bool_arr = horizons_sun_dtm_arr >= start_dtm;
    sun_long_arr = horizons_sun_long_arr(bool_arr);
    sun_dtm_arr = horizons_sun_dtm_arr(bool_arr);
    
    hour_angle = (ltim - 12)*(360/24);
    sun_long = long - hour_angle;
    sun_long = mod(sun_long, 360);
    
    sun_long_diff_arr = sun_long_arr - sun_long;
    % identify index where sign in long diff changes implying a solution exists
    interp_idx = find(sun_long_diff_arr(1:end-1) > 0 & sun_long_diff_arr(2:end) < 0);
    
    % deal with edge cases where value loops from 0 -> 360
    if sun_long > max(sun_long_arr) || sun_long < min(sun_long_arr)
        [~,interp_idx] = find(diff(sun_long_arr)>0);
    end
    
    try
        interp_idx = interp_idx(1); % want first instance after start_dtm
    catch me
        me
        interp_idx = find(sun_long_diff_arr(1:end-1) > 0 & sun_long_diff_arr(2:end) < 0)
        [~,interp_idx] = find(diff(sun_long_arr)>0)
        sum(bool_arr(:))
        
    end
    
    % deal with cases at either end of data
    if interp_idx < 2
        interp_idx = 2;
    elseif interp_idx >= numel(sun_long_arr)
        interp_idx = numel(sun_long_arr) - 1;
    end
    
    % interpolate to find more exact value around idx
    long_interp_arr = sun_long_arr(interp_idx-1:interp_idx+1);
    dtm_interp_arr = sun_dtm_arr(interp_idx-1:interp_idx+1);
    
    output_dtm(ltim_idx) = interp1(long_interp_arr, dtm_interp_arr, sun_long);
    
    % fallback
    if isnat(output_dtm(ltim_idx))
        output_dtm(ltim_idx) = sun_dtm_arr(interp_idx);
    end
end

end