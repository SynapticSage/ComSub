function plot_statistic_correlation(selected, efizz, field, statistic, varargin)
    % Parse the optional input arguments
    ip = inputParser;
    ip.addParameter('saveloc', '', @ischar);
    ip.addParameter('savetitle', '', @ischar);
    ip.parse(varargin{:});
    p = ip.Results;

    % Handle the missing data for the specific field
    selected.efizz.(field) = fillmissing(selected.efizz.(field), 'Constant', 0);
    interp_data = interp1(selected.efizz.t, selected.efizz.(field), selected.pattern.X_time);

    concatenated_matrix = [selected.pattern.u, selected.pattern.v, selected.pattern.us, selected.pattern.vs];
    
    % Calculate correlation
    [correlation_matrix,pval] = corr(interp_data, concatenated_matrix, 'Rows', 'complete');

    patterns = {'u', 'v', 'us', 'vs'};
    pattern_indices = [1, size(selected.pattern.u, 2);...
                       size(selected.pattern.u, 2)+1, size(selected.pattern.u, 2)+size(selected.pattern.v, 2);...
                       size(selected.pattern.u, 2)+size(selected.pattern.v, 2)+1, size(selected.pattern.u, 2)+size(selected.pattern.v, 2)+size(selected.pattern.us, 2);...
                       size(selected.pattern.u, 2)+size(selected.pattern.v, 2)+size(selected.pattern.us, 2)+1, size(concatenated_matrix, 2)];
    colors = [[0.7 0.7 1]; [0.4 0.4 1]; [0.9 0.5 0.5]; [1 0.2 0.2]];

    % Bootstrapping
    nBoot = 1000; % Number of bootstrap samples
    alpha = 0.05; % For 95% CI
    correlation_matrix(pval < alpha) = NaN;
    % Bootstrapping the entire correlation matrix
    nBoot = 500; % Number of bootstrap samples
    boot_corr_matrices = zeros(size(correlation_matrix, 1), size(correlation_matrix, 2), nBoot);

    for b = 1:nBoot
        % Resample indices with replacement
        resample_idx = randi(length(interp_data), length(interp_data), 1);

        % Resample data
        boot_interp_data = interp_data(resample_idx,:);
        boot_concat_matrix = concatenated_matrix(resample_idx, :);

        % Calculate correlation for the bootstrapped sample
        [c,p] = corr(boot_interp_data, boot_concat_matrix, 'Rows', 'complete');
        c(p < alpha) = nan;
        boot_corr_matrices(:,:,b) = c;
    end


    hold on;

    for i = 1:4
        % Extract data for the current pattern
        boot_data = squeeze(boot_corr_matrices(:, pattern_indices(i, 1):pattern_indices(i, 2), :));
        boot_means = feval(statistic, boot_data, 2);
        ci = prctile(boot_means, [100*alpha/2, 100*(1-alpha/2)], 2);

        f = efizz.f;
        keyboard
        y = mean(boot_means, 2);
        plot(f, y, 'Color', colors(i, :), 'LineWidth', 2);
        plots.fill_curve(f, [ci(:,1), y], colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plots.fill_curve(f, [y, ci(:,2)], colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Draw dashed line for zero-mean correlation
    plot(xlim, [0 0], '--k');
    hold off;

    % Additional plot customizations
    title(['Summary of ' func2str(statistic) ' Correlation']);
    ylabel('Correlation Value');
    xlabel('Frequency');
    legend(patterns, 'Location', 'northwest');
    grid on;

    % Save the plot if necessary
    if ~isempty(p.saveloc)
        saveas(gcf, fullfile(p.saveloc, p.savetitle + ".fig"));
        saveas(gcf, fullfile(p.saveloc, p.savetitle + ".pdf"));
        saveas(gcf, fullfile(p.saveloc, p.savetitle + ".png"));
        save(fullfile(p.saveloc, p.savetitle + ".mat"), 'correlation_matrix');
    end
end

