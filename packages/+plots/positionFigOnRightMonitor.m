function positionFigOnRightMonitor(f, figWidth, figHeight)
if nargin < 1 || isempty(f)
    f = gcf;
end
if nargin < 2 || isempty(figWidth)
    figWidth = 500;
end
if nargin < 3 || isempty(figHeight)
    figHeight = 500;
end

primaryScreenSize = get(0, 'ScreenSize'); % [left, bottom, width, height]

widthOfPrimaryMonitor = primaryScreenSize(3);

% Position the figure on the right monitor
set(f, 'Position', [widthOfPrimaryMonitor, 100, figWidth, figHeight]);
