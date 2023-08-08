function resultTable = ccatime(Components_overall, varargin)
    % Parse optional input arguments
    p = inputParser;
    addOptional(p, 'efizz', struct(), @isstruct);
    addOptional(p, 'Option', struct(), @isstruct);
    addOptional(p, 'behavior', table(), @istable);
    addParameter(p, 'behaviorColumns', {}, @iscellstr);
    addParameter(p, 'N', 3, @isnumeric); % Number of components to extract
    parse(p, varargin{:});

    efizz = p.Results.efizz;
    Option = p.Results.Option;
    behavior = p.Results.behavior;
    behaviorColumns = p.Results.behaviorColumns;
    N = p.Results.N;

    % Initialize an empty table
    resultTable = table();

    consts = option.constants;
    directions = consts.directions;

    % Loop over all components
    for i = progress(1:size(Components_overall, 1), 'Title', 'CCA time')
        for j = progress(1:size(Components_overall, 2), 'Title', 'pattern')

            % Get the cca field from the current component
            cca = Components_overall(i, j).cca;
            if isempty(cca)
                continue
            end

            % Extract the first 3 U x V components
            U = cca.u(:, 1:N);
            V = cca.v(:, 1:N);
            pattern = Option.patternNames(j);
            genH = Option.genH_name;
            direction = directions(i);
            pattern = repmat(pattern, size(U, 1), 1);
            genH = repmat(genH, size(U, 1), 1);
            direction = repmat(direction, size(U, 1), 1);

            % Create a table from the U and V data
            UVTable = array2table([U, V], 'VariableNames',...
                {'U1', 'U2', 'U3', 'V1', 'V2', 'V3'});
            strtable = array2table([pattern, genH, direction], 'VariableNames',...
                {'pattern', 'genH', 'direction'});
            UVTable = [UVTable, strtable];

            % Repeat the scalar fields to match the size of the U and V data
            scalarFields = structfun(@num2str,...
            nd.flexrmfield(Components_overall(i, j), {'X_source', 'X_target',...
            'rankRegress', 'factorAnalysis', 'jpecc', 'cca',...
            'event_anal','index_source','index_target','source','target',...
            'X_time'}), 'UniformOutput', false);
            % Create an array of scalar field values
            scalarVals        = string(struct2cell(scalarFields));
            scalarFields      = string(fieldnames(scalarFields));
            scalarVals        = num2cell(scalarVals);
            scalarFieldsTable = table(scalarVals{:}, 'VariableNames', scalarFields);
            scalarFieldsTable = repmat(scalarFieldsTable, size(U, 1), 1);

            % Convert cell columns to string columns
            for k = 1:size(scalarFieldsTable, 2)
                if iscell(scalarFieldsTable{:, k})
                    scalarFieldsTable{:, k} = string(scalarFieldsTable{:, k});
                end
            end

            % Add time column if X_time is not empty
            if ~isempty(Components_overall(i, j).X_time)
                time = Components_overall(i, j).X_time(:);
                timeTable = array2table(time, 'VariableNames', {'time'});
                % Concatenate the U V table, the time table, and the scalar fields table
                tempTable = [UVTable, timeTable, scalarFieldsTable];
            else
                % Concatenate the U V table and the scalar fields table
                tempTable = [UVTable, scalarFieldsTable];
            end

            % Add spectral properties at these times
            if ~isempty(Option)
                tempTable.animal = repmat(Option.animal, size(U, 1), 1);
                for k = 1:size(Option.frequenciesPerPattern, 1)
                    band = Option.frequenciesPerPattern(k, :);
                    idx = efizz.f >= band(1) & efizz.f <= band(2);
                    % Default spectral fields
                    fields = {'S1', 'S2', 'Cavg', 'wpli_avg'};
                    for l = 1:numel(fields)
                        disp("Adding spectral field: " + fields{l});
                        if contains(fields{l}, 'wpli')
                            efizz.(fields{l})(isnan(efizz.(fields{l}))) = 0;
                        end
                        field = fields{l};
                        data = efizz.(field)(:, idx);
                        data = mean(data, 2);
                        % Interpolate the data if it is not aligned to the times
                        if size(data, 1) ~= size(U, 1)
                            data = interp1(efizz.t, data, time, 'nearest');
                        end
                        name = field + Option.patternNames(k);
                        tempTable.(name) = data;
                    end
                end

                % Add animal name as a column
                tempTable = addvars(tempTable,...
                repmat({efizz.animalnames}, size(U, 1), 1),...
                'NewVariableNames', {'animalnames'});
            end

            % Add behavior columns
            if ~isempty(behavior) && ~isempty(behaviorColumns)
                for k = 1:numel(behaviorColumns)
                    column = behaviorColumns{k};
                    if ismember(column, behavior.Properties.VariableNames)
                        data = behavior.(column);
                        % Interpolate the data if it is not aligned to the times
                        if size(data, 1) ~= size(U, 1)
                            dataclass = class(data);
                            data = double(data);
                            data = interp1(behavior.time, data, time, 'nearest');
                            if dataclass == "logical"
                                data(isnan(data)) = false;
                                data = logical(data);
                            else
                                data = cast(data, dataclass);
                            end
                        end
                        disp("Adding behavior column: " + column + " with class " + class(data));
                        tempTable = addvars(tempTable, data, 'NewVariableNames', {column});
                    end
                end
            end

            % Append the tempTable to the resultTable
            resultTable = [resultTable; tempTable];  %#ok<AGROW>
        end
    end


if ~isempty(Option)
    writetable(resultTable,...
        figuredefine("tables", Option.animal + Option.genH_name + "ccatime.csv"),...
        'WriteRowNames', true);
else
    writetable(resultTable, figuredefine("tables", "ccatime.csv"),...
        'WriteRowNames', true);
end
