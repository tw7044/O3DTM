function colors = log_colormap(input_colors, power, value, tight_colorbar)
% LOG_COLORMAP transforms colormap into a logarithmic scale
% 
% colors = LOG_COLORMAP(...) sets colormap to logarithmic and returns array
% of colormap colors
%
% LOG_COLORMAP(input_colors) defines colormap to scale
%
% LOG_COLORMAP(input_colors, power) defines power to scale colormap by (1 =
% unscaled)
%
% LOG_COLORMAP(input_colors, power, value) defines point to scale colors
% from
% 
% LOG_COLORMAP(input_colors, power, value, tight_colorbar) defines if
% color value range should fit to data

if nargin == 0 || numel(input_colors) == 0
    input_colors = viridis(256);
end
if nargin < 2
    power = 0.5;
end
colors = NaN(256, 3);
for idx = 1:256
    frac = ((idx-1)/255)^power;
    interp_idx = 1 + frac*255;
    colors(idx,:) = interp1(1:256, input_colors, interp_idx);
end
if nargin >= 3
    lim = caxis;
    if nargin < 4 || tight_colorbar
        upper_lim = max(lim) - value;
        lower_lim = value - min(lim);
        if upper_lim > lower_lim
            upper_idx = 256;
            lower_idx = 129 - round(128*lower_lim/upper_lim);
        else
            lower_idx = 1;
            upper_idx = 128 + round(128*upper_lim/lower_lim);
        end
        colors = colors(lower_idx:upper_idx, :);
        colormap(colors)
    else
        half_lim = max(abs(value - lim));
        caxis(value + half_lim*[-1 1])
        colormap(colors)
    end
    
else
    colormap(colors)
end
end