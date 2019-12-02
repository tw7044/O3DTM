function save_report_figure(name, height, width, format)
% SAVE_REPORT_FIGURE saves figure in format used in report 
error('Report complete - don''t use this function!')

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
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
file_name = sprintf('../report/main/images/%s', name);
print(file_name, format, '-r300')
if ~strcmp(format, '-dpng')
    print(file_name, '-dpng', '-r300')
end
end