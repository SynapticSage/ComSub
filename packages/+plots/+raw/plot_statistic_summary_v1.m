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
    correlation_matrix = corr(interp_data, concatenated_matrix, 'Rows', 'complete');
    
    % Compute the desired statistic for each pattern
    num_u = size(selected.pattern.u, 2);
    num_v = size(selected.pattern.v, 2);
    num_us = size(selected.pattern.us, 2);
    num_vs = size(selected.pattern.vs, 2);

    u_stat = feval(statistic, correlation_matrix(:, 1:num_u));
    v_stat = feval(statistic, correlation_matrix(:, (num_u+1):(num_u+num_v)));
    us_stat = feval(statistic, correlation_matrix(:, (num_u+num_v+1):(num_u+num_v+num_us)));
    vs_stat = feval(statistic, correlation_matrix(:, (num_u+num_v+num_us+1):end));

    % Plotting
    hold on;

    plot(u_stat, 'Color', [0.7 0.7 1], 'LineWidth', 2);
    plot(v_stat, 'Color', [0.4 0.4 1], 'LineWidth', 2);
    plot(us_stat, 'Color', [0.9 0.5 0.5], 'LineWidth', 2);
    plot(vs_stat, 'Color', [1 0.2 0.2], 'LineWidth', 2);

    % Draw dashed line for zero-mean correlation
    plot(xlim, [0 0], '--k');
    hold off;

    % Additional plot customizations
    title(['Summary of ' func2str(statistic) ' Correlation']);
    ylabel('Correlation Value');
    xlabel('Frequency');
    legend('u', 'v', 'us', 'vs');
    grid on;

    % Save the plot if necessary
    if ~isempty(p.saveloc)
        saveas(gcf, fullfile(p.saveloc, p.savetitle + ".fig"))
        saveas(gcf, fullfile(p.saveloc, p.savetitle + ".pdf"))
        saveas(gcf, fullfile(p.saveloc, p.savetitle + ".png"))
        save(fullfile(p.saveloc, p.savetitle + ".mat"), 'u_stat', 'v_stat', 'us_stat', 'vs_stat');
    end
end

