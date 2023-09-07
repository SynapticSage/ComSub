function [avgQuantile, overallQuantiles] = test_avgQuantile(Events, fieldName)
    % Check if the field exists
    if ~isfield(Events, fieldName)
        error('Field %s does not exist.', fieldName);
    end
    
    % Extract the field
    fieldData = Events.(fieldName);
    
    % Initialize the output
    avgQuantile = zeros(1, size(fieldData, 2));
    overallQuantiles = {};
    for w = 1:length(Events.cellOfWindows)
        for col = 1:size(fieldData, 2)
            overallQuantiles{w, col} = [];
        end
    end
    
    % Loop through each column
    for col = 1:size(fieldData, 2)
        % Initialize an array to hold the quantiles
        quantiles = cell(length(Events.cellOfWindows), 1);
        
        % Sort the entire column
        sortedColumn = sort(fieldData(:, col));
        
        % Loop through each set of windows
        for w = 1:length(Events.cellOfWindows)
            for row = 1:size(Events.cellOfWindows{w},1)
                % Get the current window
                currWindow = Events.cellOfWindows{w};
                
                % Find the indices of the samples within the current window
                windowIdx = find(Events.times >= currWindow(row,1) & Events.times <= currWindow(row,2));
                
                % Extract the samples for this window
                windowSamples = fieldData(windowIdx, col);
                
                % Calculate the quantiles for these samples
                [~, windowQuantiles] = ismember(windowSamples, sortedColumn);
                windowQuantiles = windowQuantiles / length(sortedColumn);
                
                % Append to the list of quantiles
                overallQuantiles{w, col} = ...
                [overallQuantiles{w, col}; windowQuantiles];
            end
        end
        
    end
    % Calculate the average quantile for this column
    avgQuantile = cellfun(@mean, overallQuantiles);
end

