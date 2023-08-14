function [updatedTimes, ranges] = removeDataGaps(originalTimes, ranges, time_gaps, varargin)
    
    ip = inputParser;
    ip.addParameter('gap_thresh', [], @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('ploton', false, @(x) islogical(x) && isscalar(x));
    ip.parse(varargin{:});
    Opt = ip.Results;

    % epsilon = 1e-10;
    assert(issorted(ranges(:)), 'ranges must be sorted');
    % Initialize the updated times with the original times
    updatedTimes = originalTimes;
    % Calculate the total number of subplots needed: one for before loop and one for each iteration
    % Create a figure for all subplots
    if Opt.ploton
        % totalPlots = size(ranges, 2);
        figure('Name','check for gaps'); clf
        tiledlayout('flow');
        % Plot the initial state of updatedTimes
        nexttile
        plot(updatedTimes);
        title('updatedTimes: Before Adjustments');
    end
    % Iterate through the identified ranges and adjust the time stamps based on time_gaps
    % ... [rest of the code above]
    % Iterate through the identified ranges and adjust the time stamps based on time_gaps
    for i = 2:size(ranges, 2)
        % Identify indices in the current range
        idx = find(updatedTimes >= ranges(1, i));
        gap = ranges(1, i) - ranges(2, i-1); % gap between start of current and end of previous
        % Adjust the timestamps for these indices
        updatedTimes(idx) = updatedTimes(idx) - gap;
        ranges(:, i:end) = ranges(:, i:end) - gap;
        % Plot the updated state of updatedTimes for this iteration
        if Opt.ploton
            nexttile
            plot(updatedTimes);
            title(['updatedTimes: After Adjustment ', num2str(i-1)]);
        end
    end


    if ~isempty(Opt.gap_thresh)
        if (Opt.gap_thresh < 0)
            warning('thresh_gap must be positive'); end
        if all(diff(updatedTimes) <= 0)
             warning('updatedTimes must be monotonically increasing');
             d = diff(updatedTimes);
             disp("Sample violations: " + d(d<0));
        end
        if all(diff(updatedTimes) >= Opt.gap_thresh)
            warning('thresh_gap must be greater than the largest gap in updatedTimes');
        end
    end

end

