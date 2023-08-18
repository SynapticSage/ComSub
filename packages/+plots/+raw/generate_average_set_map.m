function generate_average_set_heatmap(selected, selected_rows, sets_wanna_plot, varargin)
% generate_average_set_heatmap: Generates a heatmap of the average
    
% Input parsing
p = inputParser;
addRequired(p, 'selected');
addRequired(p, 'selected_rows');
addRequired(p, 'sets_wanna_plot');
addParameter(p, 'split_by_reward', false);
addParameter(p, 'grid_res', 100);
addParameter(p, 'sgtitlePrepend', '');
addParameter(p, 'useRollingStd', false); % new option for using rolling std
addParameter(p, 'rollingStdWindow', 13); % default window size for rolling std
addParameter(p, 'saveloc', '');
addParameter(p, 'savetitle', '');
parse(p, selected, selected_rows, sets_wanna_plot, varargin{:});

split_by_reward = p.Results.split_by_reward;
grid_res = p.Results.grid_res;

% Grid discretization using behavior data from selected_rows
x_indices_all = round(rescale(selected_rows.X, 1, grid_res));
y_indices_all = round(rescale(selected_rows.Y, 1, grid_res));

% Determine the number of subplots
disp("Computing size")
numSets = length(sets_wanna_plot);
subplot_dim1 = ceil(sqrt(numSets));
subplot_dim2 = floor(sqrt(numSets));
if split_by_reward
    subplot_dim2 = subplot_dim2 * 2; % double the columns
    while subplot_dim1 * subplot_dim2 < numSets*2
        disp("Increasing subplot_dim1")
        subplot_dim1 = subplot_dim1 + 1;
    end
else
    while subplot_dim1 * subplot_dim2 < numSets
        disp("Increasing subplot_dim1")
        subplot_dim1 = subplot_dim1 + 1;
    end
end

% Setup the figure
sgtitle_text = 'Trajectory Heatmaps';
if ~isempty(p.Results.sgtitlePrepend)
    sgtitle_text = char(string(p.Results.sgtitlePrepend) + newline + sgtitle_text);
end
f = fig(sgtitle_text);
disp("Plotting")
for set_idx = 1:numSets
    % Extract the current set_to_plot parameters
    set_to_plot = sets_wanna_plot{set_idx};
    time_name   = set_to_plot(3);
    data_name   = set_to_plot(1);
    comp_name   = double(set_to_plot(2));
    set_to_plot_title = data_name + ", " + comp_name + ' vs ' + time_name;
    set_to_plot_saveappend = data_name + "_" + comp_name + '_vs_' + time_name;
    % Compute rolling standard deviation if the flag is set
    if p.Results.useRollingStd
        selected.pattern.(data_name) = movstd(selected.pattern.(data_name), p.Results.rollingStdWindow, 'omitnan');
    end
    data_for_coloring = interp1(selected.pattern.(time_name), selected.pattern.(data_name)(:,comp_name), selected.behavior.time, 'nearest','extrap');
    extrap_locations = isnan(interp1(selected.pattern.(time_name), selected.pattern.(data_name)(:,comp_name), selected.behavior.time, 'nearest'));
    data_for_coloring(extrap_locations) = NaN;
    
    % Form title string
    title_str = string(data_name) + ", " + comp_name + ' vs ' + time_name;

    % Check which trajectories to plot based on the set_to_plot conditions
    inds_set = true(height(selected_rows), 1); % This might need further customization based on conditions

    % Create the heatmap grid based on rewarded value if split_by_reward is true
    if split_by_reward
        rewarded_values = unique(selected_rows.rewarded);
        for r = 1:length(rewarded_values)
            data_grid = zeros(grid_res);
            heatmap_grid = zeros(grid_res);
            inds = inds_set & selected_rows.rewarded == rewarded_values(r);
            for i = find(inds)'
                if isnan(data_for_coloring(i))
                    continue;
                end
                data_grid(y_indices_all(i), x_indices_all(i)) = data_grid(y_indices_all(i), x_indices_all(i)) + data_for_coloring(i);
                heatmap_grid(y_indices_all(i), x_indices_all(i)) = heatmap_grid(y_indices_all(i), x_indices_all(i)) + 1;
            end

            % Normalize and plot
            heatmap_grid = heatmap_grid / sum(heatmap_grid(:));
            heatmap_grid = data_grid ./ heatmap_grid;
            subplot(subplot_dim1, subplot_dim2, (set_idx-1)*2 + r); 
            imagesc(heatmap_grid);
            if ~p.Results.useRollingStd
                cm = [cmocean('balance',1024); 0 0 0];
                colormap(cm);
                clim([-nanmax(abs(heatmap_grid(:))), nanmax(abs(heatmap_grid(:)))]);
            else
                cm = [cmocean('thermal',1024); 0 0 0];
                colormap(cm);
            end
            title([title_str, ' Rewarded:', num2str(rewarded_values(r))]);
            colorbar;
            axis equal;
            axis tight;
        end
    else
        heatmap_grid = zeros(grid_res);
        data_grid    = zeros(grid_res);
        inds = inds_set;
        for i = find(inds)'
            if isnan(data_for_coloring(i))
                continue;
            end
            data_grid(y_indices_all(i), x_indices_all(i)) = data_grid(y_indices_all(i), x_indices_all(i)) + data_for_coloring(i);
            heatmap_grid(y_indices_all(i), x_indices_all(i)) = heatmap_grid(y_indices_all(i), x_indices_all(i)) + 1;
        end

        % Normalize and plot
        heatmap_grid = heatmap_grid / sum(heatmap_grid(:));
        heatmap_grid = data_grid ./ heatmap_grid;
        subplot(subplot_dim1, subplot_dim2, set_idx);
        imagesc(heatmap_grid);
        if ~p.Results.useRollingStd
            cm = [cmocean('balance',1024); 0 0 0];
            colormap(cm);
            clim([-nanmax(abs(heatmap_grid(:))), nanmax(abs(heatmap_grid(:)))]);
        else
            cm = [cmocean('thermal',1024); 0 0 0];
            colormap(cm);
        end
        title(title_str);
        colorbar;
        axis equal;
        axis tight;
    end
    if ~isempty(p.Results.saveloc)
        disp("Saving " + fullfile(p.Results.saveloc, p.Results.savetitle + set_to_plot_saveappend + ".mat"))
        save(fullfile(p.Results.saveloc, p.Results.savetitle + set_to_plot_saveappend + ".mat"), 'heatmap_grid', 'data_grid', 'title_str', 'set_to_plot');
    end
end
set(gcf, 'Position', get(0, 'Screensize'));
if ~isempty(p.Results.savetitle)
    disp("Saving " + fullfile(p.Results.saveloc, p.Results.savetitle + ".fig"))
    saveas(gcf, fullfile(p.Results.saveloc, p.Results.savetitle + ".fig"));
    saveas(gcf, fullfile(p.Results.saveloc, p.Results.savetitle + ".pdf"));
    saveas(gcf, fullfile(p.Results.saveloc, p.Results.savetitle + ".png"));
end
end

