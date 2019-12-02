function colors = step_colormap(input_colors, power, value, tight_colorbar, ax_h)
% STEP_COLORMAP transforms colormap into a step-like function
%
% colors = STEP_COLORMAP(...) sets colormap to step colormap and returns
% matrix of colors
%
% STEP_COLORMAP(input_colors, power) specifies colormap to transform and
% power to define step steepness (smaller = steeper)
%
% STEP_COLORMAP(input_colors, power, value) defines the central value of
% the colormap (i.e. where the step is)
%
% STEP_COLORMAP(input_colors, power, value, tight_colorbar) defines if
% colorbar should show full range of colormap values or fit data
%
% STEP_COLORMAP(input_colors, power, value, tight_colorbar, ax_h) defines
% axis handle to apply colormap to

if nargin == 0 || numel(input_colors) == 0
    input_colors = viridis(256);
end
if nargin < 2
    power = 0.5;
end
colors = NaN(256, 3);
for idx = 1:256
    frac = (abs(idx-128)/128)^power*sign(idx-128);
    interp_idx = 128+frac*127;
    colors(idx,:) = interp1(1:256, input_colors, interp_idx);
end
if nargin >= 3
    if nargin >= 5
        lim = caxis(ax_h);
    else
        lim = caxis;
    end
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
    else
        half_lim = max(abs(value - lim));
        caxis(value + half_lim*[-1 1])
    end
end
if nargin >= 5
    colormap(ax_h, colors)
else
    colormap(colors)
end
end