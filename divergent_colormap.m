function divergent_colormap(limit_value, reversed, centre_value)
% DIVERGENT_COLORMAP changes colormap to colormap for differences with blue
% for negative values and red for positive values (with color intensity
% proportional to the difference magnitude).
%
% DIVERGENT_COLORMAP sets colormap to red-yellow-blue with limits taken
% from caxis
%
% DIVERGENT_COLORMAP(limit_value) sets caxis limits as [-limit_value,
% +limit_value]
% 
% DIVERGENT_COLORMAP(limit_value, reversed) flips order of colormap
%
% DIVERGENT_COLORMAP(limit_value, reversed, centre_value) sets caxis limits
% to centre_value + [-limit_value, +limit_value]


if nargin < 2
    reversed = false;
end
if reversed
    colormap(flipud(spectral))
else
    colormap(spectral)
end
if nargin == 0
    caxis('auto')
    lim = max(abs(caxis));
end
if nargin > 0 && numel(limit_value > 0)
    lim = limit_value;
end
caxis([-lim, +lim])
if nargin == 3
    caxis(centre_value + limit_value*[-1 1])
end
end