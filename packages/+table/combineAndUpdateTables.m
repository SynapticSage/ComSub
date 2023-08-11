function out = combineAndUpdateTables(regex, finalFile)
% COMBINEANDUPDATETABLES Combine and update tables
%   COMBINEANDUPDATETABLES(regex, finalFile) combines and updates tables
%   that match the regular expression regex and saves the combined table to
%   finalFile. The tables must have a hash column.
%
%   Example:
%       combineAndUpdateTables('*.mat', 'finalTable.mat')

    % Load the existing table if it exists
    if ~endsWith(finalFile, '.mat')
        finalFile = finalFile + ".mat";
    end
    if exist(finalFile, 'file')
        out = load(finalFile);
        fields = fieldnames(out);
        assert(numel(fields) == 1);
        out = out.(fields{1});
    else
        out = table();
    end
    out = table.hackyconversions(out);

    % Get a list of the new table files
    newFiles = dir(codedefine("DATA_TABLES", regex));
    if ~contains(regex, '*')
        warning('No wildcards found in regex');
    end
    if isempty(newFiles)
        warning('No new table files found');
        return
    end

    % Loop over each new table file
    for i = 1:length(newFiles)

        % Load the new table
        newTable = load(fullfile(newFiles(i).folder, newFiles(i).name));
        fields = fieldnames(newTable);
        assert(numel(fields) == 1);
        newTable = newTable.(fields{1});
        newTable = table.hackyconversions(newTable);

        % Align the existing table and the new table
        [out, newTable] = table.alignTables(out, newTable);
        
        % Loop over each row in the new table
        for j = 1:height(newTable)
            % Find if the hash value exists in the existing table
            idx = find(strcmp(out.hash, newTable.hash(j)));
            
            if isempty(idx)
                % If the hash value does not exist, append the new row
                % out = [out; newTable(j, :)];
                out = table.flexvertcat(out, newTable(j, :));
            else
                % If the hash value exists, replace the existing row
                for f = setdiff(fieldnames(newTable)', {'Properties', 'Row', 'Variables'})
                    try
                        out.(f{1})(idx) = newTable.(f{1})(j);
                    catch
                        warning('Could not update %s', f{1});
                        keyboard
                    end
                end
            end
        end
    end

    disp("Combined " + num2str(length(newFiles)) + " tables into " + finalFile);
    disp("Final table")
    disp(out)
    
    % Save the updated table
    finalFile = replace(finalFile, '.mat', '');
    out = struct(finalFile, out);
    save(finalFile, '-struct', 'out');
    
    % Delete all the tables used
    for i = 1:length(newFiles)
        delete(fullfile(newFiles(i).folder, newFiles(i).name));
    end
end
