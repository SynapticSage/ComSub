function plot_event_values(out_struct, varargin)
%PLOT_EVENT_VALUES Plot the event values
% see cca.event_analysis for details on the out_struct

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
ip.addParameter('sample', 800, @(x) isnumeric(x) && isscalar(x));
ip.parse(varargin{:});
Opt = ip.Results;

% Get the event_u_values and event_v_values
event_u_values = out_struct.event_u_values;
event_v_values = out_struct.event_v_values;

% Get the sizes of the 3D matrices
[numRows, numSamples, numCols] = size(event_u_values);

% Create a new figure
figure;

% Create a 3D tiled layout
Nh = numRows;
Nw = numCols;
if Nh == 0 || Nw == 0
    disp("No event values to plot");
    return;
end
gap    = [0.01 0.03];
marg_h = [0.1 0.01];
marg_w = [0.01 0.01];
clear ha
ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w);

points_u = cell(numRows, numCols);
points_v = cell(numRows, numCols);
% Loop over the rows and columns
for row = 1:numRows
    for col = 1:numCols
        % Select the subplot for this row and column
        t = ha((row - 1) * numCols + col);
        % Set axes
        axes(t);

        % Get the u and v values for this row and column
        u_values = squeeze(event_u_values(row, :, col));
        v_values = squeeze(event_v_values(row, :, col));
        inds = randperm(numSamples);
        u_values = u_values(inds);
        v_values = v_values(inds);
        points_u{row, col} = u_values;
        points_v{row, col} = v_values;

        % Plot the u and v values
        scatter(u_values, v_values);
        axis equal;
        % hline at 0
        xline(0);
        % vline at 0
        yline(0);
        alpha(0.25);

        % X=Y line with a text label "communication subspace" along
        % the line (45 degrees, 1.5 units from the origin)
        line([-3 3], [-3 3], 'Color', 'red', 'LineStyle', '--');
        tt=text(-1, -1.3, 'communication subspace', 'Color', 'red');
        tt.Rotation = 45;


        % Set the title
        title(sprintf('Pattern: %d, Component: %d', row, col));
    end
end

disp("setting axes");
uq = cat(1, points_u{:});
vq = cat(1, points_v{:});
qlim_x = quantile(uq(:), [0.001 0.999]);
qlim_y = quantile(vq(:), [0.001 0.999]);

linkaxes(ha, 'xy');
set(findobj(gcf,'type','axes'),'XLim',qlim_x,'YLim',qlim_y);
set(gcf, 'Position', get(0, 'Screensize'));
pos=get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) pos(3)*1.5 pos(4)*1.5]);
if ~exist(figuredefine("cca_eventanal"), 'dir')
    mkdir(figuredefine("cca_eventanal"));
end
disp("saving");
saveas(gcf, figuredefine("cca_eventanal", sprintf('event_values%s.png', Opt.figAppend)));
saveas(gcf, figuredefine("cca_eventanal", sprintf('event_values%s.pdf', Opt.figAppend)));
saveas(gcf, figuredefine("cca_eventanal", sprintf('event_values%s.fig', Opt.figAppend)));
end
