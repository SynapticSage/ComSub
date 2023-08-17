function generate_average_set_heatmap(selected, selected_rows, sets_wanna_plot, varargin)
    
    % Input parsing
    p = inputParser;
    addRequired(p, 'selected');
    addRequired(p, 'selected_rows');
    addRequired(p, 'sets_wanna_plot');
    addParameter(p, 'split_by_reward', false);
    addParameter(p, 'grid_res', 100);
    parse(p, selected, selected_rows, sets_wanna_plot, varargin{:});
    
    split_by_reward = p.Results.split_by_reward;
    grid_res = p.Results.grid_res;

    % Grid discretization using behavior data from selected_rows
    x_indices_all = round(rescale(selected_rows.X, 1, grid_res));
    y_indices_all = round(rescale(selected_rows.Y, 1, grid_res));

    % Determine the number of subplots
    numSets = length(sets_wanna_plot);
    subplot_dim1 = ceil(sqrt(numSets));
    subplot_dim2 = floor(sqrt(numSets));
    if split_by_reward
        subplot_dim2 = subplot_dim2 * 2; % double the columns
    end

    % Setup the figure
    f = fig('Trajectory Heatmaps');

    for set_idx = 1:numSets
        % Extract the current set_to_plot parameters
        set_to_plot = sets_wanna_plot{set_idx};
        time_name = set_to_plot{1}(3);
        data_name = set_to_plot{1}(1);
        comp_name = double(set_to_plot{1}(2));
        
        % Form title string
        title_str = strcat(string(data_name) + ", " + comp_name, ' vs ', time_name);

        % Check which trajectories to plot based on the set_to_plot conditions
        inds_set = true(height(selected_rows), 1); % This might need further customization based on conditions

        % Create the heatmap grid based on rewarded value if split_by_reward is true
        if split_by_reward
            rewarded_values = unique(selected_rows.rewarded);
            for r = 1:length(rewarded_values)
                heatmap_grid = zeros(grid_res);
                inds = inds_set & selected_rows.rewarded == rewarded_values(r);
                for i = find(inds)'
                    heatmap_grid(y_indices_all(i), x_indices_all(i)) = heatmap_grid(y_indices_all(i), x_indices_all(i)) + 1;
                end

                % Normalize and plot
                heatmap_grid = heatmap_grid / max(heatmap_grid(:));
                subplot(subplot_dim1, subplot_dim2, (set_idx-1)*2 + r);
                imagesc(heatmap_grid);
                colormap(cmocean('balance'));
                title([title_str, ' Rewarded:', num2str(rewarded_values(r))]);
                colorbar;
                axis equal;
                axis tight;
            end
        else
            heatmap_grid = zeros(grid_res);
            inds = inds_set;
            for i = find(inds)'
                heatmap_grid(y_indices_all(i), x_indices_all(i)) = heatmap_grid(y_indices_all(i), x_indices_all(i)) + 1;
            end

            % Normalize and plot
            heatmap_grid = heatmap_grid / max(heatmap_grid(:));
            subplot(subplot_dim1, subplot_dim2, set_idx);
            imagesc(heatmap_grid);
            colormap(cmocean('balance'));
            title(title_str);
            colorbar;
            axis equal;
            axis tight;
        end
    end
end

