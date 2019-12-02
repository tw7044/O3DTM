function save_paper_figure(name, height, width, format)
% SAVE_PAPER_FIGURE saves figure in format used in report
%
% SAVE_PAPER_FIGURE(name) specifies target file name
%
% SAVE_PAPER_FIGURE(name, height) specifies image height (inches)
%
% SAVE_PAPER_FIGURE(name, height, width) specifies image width (inches)
%
% SAVE_PAPER_FIGURE(name, height, width, format) specifies image format

if nargin < 2 || isempty(height)
    height = 3;
end
if nargin < 3 || isempty(width)
    width = 4;
end
if nargin < 4
    format = '-dpng';
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 width height];
fig.Units = 'inches';
fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
drawnow
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
file_name = sprintf('../paper/images/%s', name);
print(file_name, format, '-r300')
end