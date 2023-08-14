function [table1, table2] = alignTables(table1, table2, varargin)
% ALIGNTABLES Aligns two tables by adding missing columns to each table
%  [table1, table2] = ALIGNTABLES(table1, table2) adds missing columns to
%  each table so that they have the same columns. The missing columns are
%  added to each table and filled with NaNs.
    ip = inputParser;
    ip.addParameter('aggressive', false, @islogical);
    ip.parse(varargin{:});
    Opt = ip.Results;

    % Get the column names of each table
    cols1 = table1.Properties.VariableNames;
    cols2 = table2.Properties.VariableNames;
    
    % Find the columns that are missing from each table
    colsToAdd1 = setdiff(cols2, cols1);
    colsToAdd2 = setdiff(cols1, cols2);
    
    % Add the missing columns to each table
    for i = 1:length(colsToAdd1)
        table1.(colsToAdd1{i}) = NaN(height(table1), 1);  % Replace NaN with a suitable default value
    end
    for i = 1:length(colsToAdd2)
        table2.(colsToAdd2{i}) = NaN(height(table2), 1);  % Replace NaN with a suitable default value
    end

    % If the hash value exists, replace the existing row
    common_columns = intersect(table1.Properties.VariableNames, table2.Properties.VariableNames);
    for k = 1:numel(common_columns)
        columnName = common_columns{k};
        if ~isa(table1.(columnName), class(table2.(columnName)))
            try
                table1.(columnName) = feval(class(table2.(columnName)), table1.(columnName));
            catch ME
                if isa(table1.(columnName), 'cell') && all(cellfun(@isempty, table1.(columnName)))
                    table1.(columnName) = repmat({feval(class(table2.(columnName)), NaN)}, height(table1), 1);
                else
                    warning('Failed to convert column %s from %s to %s: %s', ...
                            columnName, class(table1.(columnName)), class(table2.(columnName)), ME.message);
                    if Opt.aggressive
                        table1.(columnName) = [];
                        table2.(columnName) = [];
                    end
                end
            end
        end
    end

assert(isequal(sort(fieldnames(table1)), sort(fieldnames(table2))), "Non-matching columns: " + union(setdiff(fieldnames(table1), fieldnames(table2)), setdiff(fieldnames(table2), fieldnames(table1))))

end


