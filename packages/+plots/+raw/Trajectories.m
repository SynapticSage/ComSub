% Extract unique trajectories
unique_traj = unique(selected_rows.trajall);
% Randomly sample 10,000 points from the main behavior table
rand_indices = randperm(height(behavior), min(10000, height(behavior)));
sampled_X = behavior.X(rand_indices);
sampled_Y = behavior.Y(rand_indices);
xlims = quantile(sampled_X, [0.01, 0.99]);
ylims = quantile(sampled_Y, [0.01, 0.99]);
pick = 100;
disp("Picking " + pick + " trajectories to plot");
L = length(unique_traj);
randset = sort(randperm(L, min(L,pick)));

% Parameters for the behavior plot
for set_to_plot = sets_wanna_plot(:)' %FOR_SET_TO_PLOT

    % time_name = "X_time";
    % data_name = "v";
    % data_col  = 2;
    time_name = set_to_plot{1}(3);
    data_name = set_to_plot{1}(1);
    data_col  = double(set_to_plot{1}(2));
    normalized = "minmaxabs";
    background_color = [0 0 0];
    % background_color = [1 1 1];

    % Create figure and tiled layout
    f = fig(animal + " Behavior Trajectories " +  epoch_option + " " + sortprop + " trajbound " + trajbound_option + newline + " dataname: " + data_name + ", col: " + data_col);clf;
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

    % Color data based on some field in the selected structure
    data_for_coloring = interp1(selected.pattern.(time_name), selected.pattern.(data_name)(:,data_col), selected.behavior.time, 'nearest','extrap');
    extrap_locations = isnan(interp1(selected.pattern.(time_name), selected.pattern.(data_name)(:,data_col), selected.behavior.time, 'nearest'));
    % Colormap for the data (can be changed to any other colormap)
    cm=cmocean('balance', 1024);
    % Normalize the interpolated data to be within [0, 1]
    if normalized == "minmax"
        normalized_data = (data_for_coloring - min(data_for_coloring)) / (max(data_for_coloring) - min(data_for_coloring));
        clims = [min(data_for_coloring), max(data_for_coloring)];
    elseif normalized == "minmaxabs"
        normalized_data = ((data_for_coloring) / (max(abs(data_for_coloring))))/2 + 0.5
        clims = [-max(abs(data_for_coloring)), max(abs(data_for_coloring))];
    end
    % Map the normalized data to colormap indices
    color_indices = round(normalized_data * (size(cm, 1) - 1) + 1);
    % Convert indices to RGB values
    colors_for_plotting = cm(color_indices, :);
    colors_for_plotting(extrap_locations, :) = repmat([0 0 0], sum(extrap_locations), 1);

    % Loop through each unique trajectory and plot X and Y positions
    for traj_num = progress(unique_traj(randset)','Title', 'Plotting Trajectories')
        % Extract the current trajectory data
        inds_traj = selected_rows.trajall == traj_num;
        curr_traj = selected_rows(inds_traj, :);
        if ~isempty(data_for_coloring) && any(extrap_locations(inds_traj))
            % If there are any extrapolated locations, plot them in black
            warning('There are extrapolated locations in trajectory %d', traj_num);
            continue;
        end
        
        % Next tile for the current trajectory
        ax = nexttile;
        
        % Plot the sampled background trajectory
        plot(ax, sampled_X, sampled_Y, '.', 'Color', [0.5 0.5 0.5]);
        hold on;
        
        if ~isempty(data_for_coloring)
            % Plot the trajectory colored by the data
            scatter(ax, curr_traj.X, curr_traj.Y, 40, colors_for_plotting(inds_traj,:), 'filled');
            colormap(ax, cm);
            h=colorbar;
            h.Label.String = 'u';
            h.TickLabels = cellstr(string(linspace(clims(1), clims(2), 5)'));
            h.Ticks = linspace(0, 1, 5);
        else
            % Plot the trajectory without coloring
            plot(ax, curr_traj.X, curr_traj.Y, 'k');
        end
        
        % Label for the start and end of trajectory
        plot(ax, curr_traj.X(1), curr_traj.Y(1), 'go');
        plot(ax, curr_traj.X(end), curr_traj.Y(end), 'ro');
        
        % Set title
        if curr_traj.rewarded(1) == 1
            rew_str = 'Rewarded';
        else
            rew_str = 'Not Rewarded';
        end
        if curr_traj.leftright(1)
            lr_str = 'Left';
        else
            lr_str = 'Right';
        end
        title(['Epoch:' num2str(epoch_option) ' Traj: ' num2str(traj_num) newline 'TrajBound: ', num2str(curr_traj.trajbound(1)), ', ', rew_str, ', ', lr_str]);
        
        % Additional axis properties (if needed)
        xlabel('X Position');
        ylabel('Y Position');
        set(ax, 'Color', background_color);
        axis equal; % Make the X and Y axis scales the same
        ylim(ylims);
        xlim(xlims);
        grid on;
    end
    sgtitle("Trajectories " +  animal + epoch_option + " trajbound=" + trajbound_option + newline + " dataname: " + data_name + ", col: " + data_col);
    sgt = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));  % Maximize the figure window
    r = @(x) "traj_" + string(replace(replace(replace(x, " ", "_"), ":", "_"),newline,"_"));
    saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.png'));
    saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.pdf'));
    saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.fig'));

end % FOR_SET_TO_PLOT
