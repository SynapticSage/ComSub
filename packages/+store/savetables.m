function [RunsSummary, DetailedRunsSummary] = savetables(Events, Patterns, Option)

% Ready the path
table_folder = fullfile(codedefine(), 'DATA_TABLES');
if ~exist(table_folder, 'dir')
    mkdir(table_folder);
end
if ~ismember(table_folder, path)
    addpath(table_folder);
end

%% RunsSummary: Abbreviated summary of results for runs
tableFolder = fullfile(string(codedefine()), "DATA_TABLES");
tabname     = @(x) fullfile(tableFolder, x);
if exist(tabname("RunsSummary" + Option.tableAppend + ".mat"), 'file') 
    load(tabname("RunsSummary" + Option.tableAppend + ".mat"), 'RunsSummary');
else
    RunsSummary = table();
end

%% DetailedRunsSummary: Summary of information for runs
if exist(tabname("DetailedRunsSummary" + Option.tableAppend + ".mat"), 'file') 
     load(tabname("DetailedRunsSummary" + Option.tableAppend + ".mat"), 'DetailedRunsSummary');
else
     warning("no existing DetailedRunsSummary table")
     DetailedRunsSummary = table();
end

RunsSummary = table.hackyconversions(RunsSummary);
DetailedRunsSummary = table.hackyconversions(DetailedRunsSummary);

% Identifying information about this options set and date of run
hash = store.gethash(Option);
timestamp = string(datetime());

% Determine information to add to table
Optim=params.getOptimizationParams(Patterns,Events,Option);
Optimtable=struct2table(Optim);
Optimtable.timestamp     = timestamp;
Optimtable.hash          = hash;
Optimtable.numWindowsCut = Events.nWindows;
Optimtable.cutoffs       = Events.cutoffs;

% Create a table of options and table of patterns
Optiontable  = struct2table(Option, 'AsArray', true);
Patterntable = query.getPatternTable(Patterns, Option);
Optiontable = table.hackyconversions(Optiontable);
Patterntable = table.hackyconversions(Patterntable);

% Combine option columns with hash and date
tablerow = [Optiontable, Optimtable];
% Combine those with all rows of the pattern table
tablecontent = util.table.flexibleColumnCat(Patterntable, ...
		repmat(tablerow, height(Patterntable), 1));

%% Check and Hash
if ~isempty(RunsSummary)  &&any(contains(RunsSummary.hash, hash))
     RunsSummary(contains(RunsSummary.hash, hash), :) = []; % Delete any rows that contain the current hash
     DetailedRunsSummary(contains(DetailedRunsSummary.hash, hash), :) = [];
     disp("already computed before, rehashing to the same location");
% New options:    Append row
else
     disp("new results  --not in existing table")
end

% --------------------------------------
% Append the new results and posrpocess
% --------------------------------------
if istable(DetailedRunsSummary)
     old_height = height(DetailedRunsSummary);
else
     old_height = 0;
end
fields_equal = isequal(fields(DetailedRunsSummary), fields(tablecontent))
if  old_height ~= 0 ||~fields_equal
     % ~isequal(fields(DetailedRunsSummary), fields(tablecontent))
     DetailedRunsSummary = table.addNewColumn(DetailedRunsSummary, tablecontent);
     RunsSummary         = table.addNewColumn(RunsSummary, tablerow);
else
     DetailedRunsSummary = [DetailedRunsSummary; tablecontent];
     RunsSummary         = [RunsSummary; tablerow];
end
assert(height(DetailedRunsSummary) > old_height, "appending failed!");
DetailedRunsSummary = table.postprocessSummaryTable(DetailedRunsSummary);
RunsSummary         = table.postprocessSummaryTable(RunsSummary);


%% ------------- Save ----------------------------
% save the tables
save(tabname("RunsSummary" + Option.tableAppend),         "RunsSummary",         '-v7.3');
save(tabname("DetailedRunsSummary" + Option.tableAppend), "DetailedRunsSummary", '-v7.3');
