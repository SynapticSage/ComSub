function fixRunsSummary(tableFile, columnName, optionAttribute, varargin)
    % Input parsing
    p = inputParser;
    addRequired(p, 'tableFile', @ischar);
    addRequired(p, 'columnName', @ischar);
    addRequired(p, 'optionAttribute', @ischar);
    addParameter(p, 'startDate', '01-Jan-1970', @ischar);
    addParameter(p, 'defaultValue', [])
    parse(p, tableFile, columnName, optionAttribute, varargin{:});
    
    startDate = datetime(p.Results.startDate, 'InputFormat', 'dd-MMM-yyyy');
    defaultValue = p.Results.defaultValue;
    
    % Load RunsSummary table from the file
    originalpath = which(tableFile);
    load(tableFile);

    % Convert timestamp column to datetime
    oldts=RunsSummary.timestamp;
    RunsSummary.timestamp = datetime(RunsSummary.timestamp, 'InputFormat', 'dd-MMM-yyyy');
    for i = 1:length(oldts)
        RunsSummary.timestamp(i) = datetime(oldts(i));
    end

    % Filter rows based on startDate
    rowsToUpdate = RunsSummary.timestamp > startDate;
    filteredSummary = RunsSummary(rowsToUpdate, :);
    
    % Iterate through each hash in the filtered RunsSummary
    nRows = height(filteredSummary);
    first = true;
    for i = progress(1:nRows, 'Title', 'Fixing RunsSummary')
        hashValue = filteredSummary.hash(i);
        matFileName = char(strjoin([hashValue, '.mat'],''));

        % Load the corresponding MAT file if it exists
        if exist(matFileName, 'file')
            loadedData = matfile(matFileName);
            rowToUpdate = RunsSummary.hash == hashValue;

            % if first
            %     keyboard
            % end

            % Check if Option struct exists in the MAT file
            if ismember('Option', who(loadedData))
                Option = loadedData.Option;

                % Check if the attribute exists in Option
                if isfield(Option, optionAttribute)
                    disp("Fixing value to " + Option.(optionAttribute) + " for " + hashValue);
                    RunsSummary.(columnName)(rowToUpdate) = Option.(optionAttribute);
                else
                    warning(['Option attribute "', optionAttribute, '" not found in ', matFileName]);
                    if ~isempty(defaultValue)
                        disp("Fixing value to " + defaultValue + " for " + hashValue);
                        RunsSummary.(columnName)(rowToUpdate) = defaultValue;
                    end
                end
            else
                warning(['Option struct not found in ', matFileName]);
            end
            first = false;
        else
            warning([matFileName, ' does not exist.']);
        end
    end

    % Save the updated RunsSummary back to the file
    disp("Check the updated RunsSummary table in the workspace.");
    keyboard
    save(tableFile)
end

