% paths
tic
if ~exist('csubpath', 'var')
    commsubspaceToPath;
    addpath(hashdefine())
    csubpath=true;
end
const = option.constants();
%% Combine tables if they exist
table.combineAndUpdateTables("RunsSummary_*", "RunsSummary");
load("RunsSummary.mat", "RunsSummary");

% show all genH_name == wpli columns
RunsSummary(contains(RunsSummary.genH_name, "wpli"),... 
["animal", "genH_name","preProcess_zscore", "hash", "timestamp"])

%% load
multi_epoch = false; % usually first 3, last 3 epochs
zscr        = true; % zscored or not
midpattern  = true; % midpattern or not
load("RunsSummary.mat");
disp(" ---->  Multi epoch: " + multi_epoch)
disp(" ---->  Zscored: " + zscr)
if zscr == false && midpattern == false
    figuredefine("-permfolder", "zscore=false");
    figuredefine("-creation", true);
elseif midpattern == true && zscr == false
    figuredefine("-permfolder", "zscore=false_midpattern=true");
    figuredefine("-creation", true);
elseif midpattern == true && zscr == true
    figuredefine("-permfolder", "midpattern=true");
    figuredefine("-creation", true);
elseif midpattern == false && zscr == true
    % nothing, default state
else
    figuredefine("-clearpermfolder")
end
if zscr
    disp("USING ZSCORE");
    disp("USING ZSCORE");
    disp("USING ZSCORE");
end
if midpattern
    disp("USING MIDPATTERN");
    disp("USING MIDPATTERN");
    disp("USING MIDPATTERN");
end
pause(0.25);

% ----------------------
% TABLE : ACQUIRE RUNS
% ----------------------
% Determine keys to use : you can use this string to arbitrarily select rows
%   each item of the filtstring is a property to select. $x pulls the x
%   column and applies the test shown
filtstring = ...
["ismember($animal, [""JS21"",""ZT2"",""ER1"",""JS14"",""JS13"",""JS17"",""JS15""])",... % RY added 2023
       ..."$spikeBinSize==0.15",...
       "$preProcess_zscore=="+zscr,...
       ..."$numPartition==50",...
       "$quantileToMakeWindows == 0.85",...
       "$midpattern=="+midpattern,...
       ..."arrayfun(@(x)isequal($winSize(x,:), [0,0.3]), 1:size($winSize,1))'" ... 
]; % WARNING: reason why JS15 missing?
% Get the proper keys
matching_runs = query.getHashed_stringFilt(RunsSummary, filtstring,...
                'mostrecent', ["animal", "generateH","preProcess_zscore"]);
disp("Number of matches found: " + height(matching_runs))
matching_runs(:, ["animal", "generateH", "preProcess_zscore", "hash", "timestamp"]) 

%%

if multi_epoch
    keys = [];
    for i = 1:height(matching_runs)
        if numel(matching_runs(i,:).epochsSelected{1} ) == 3 % change this to 3 for 3 early and 3 late
            keys = [keys, matching_runs(i,:).hash];
        end
    end
else
    keys = matching_runs.hash;
end

%% Load keyed save files from local or server
% Load up combined pattern data
localmode = true; % set to false to load from server
if localmode
    % Load locally
    Out = query.combinePatterns(keys, ...
    'RunsSummary', RunsSummary, ...
     'verbose', ["animal", "generateH", "preProcess_zscore"]);
    Patterns = Out.Patterns;
    Option   = Out.Option;
    clear Out
else
    % Grab from server and load
    [Patterns, otherData] = query.pattern.remote_LoadAndCombine(keys, 'verbose', ...
                                ["animal", "generateH", "preProcess_zscore"]);
    Option = otherData{1}.Option;
end

assert(~isempty(Option), "Data is empty -- downstream code will fail")
Patterns = nd.merge(Patterns, Option, 'only', {'animal', 'generateH', 'genH_name'}, ...
                                            'broadcastLike', true,...
                                             'overwrite', true);
[nDataset, nPartition, ~, nResult] = size(Patterns);
nPatterns = nResult/2;
szPatterns = size(Patterns);
Patterns = squeeze(Patterns);
Patterns = nd.dimLabel(Patterns, 1:ndims(Patterns), ...
            ["iDataset", "iPartition", "iDirection", "iPatterns"]);
Patterns = reshape(Patterns, szPatterns);
Patterns = util.type.castefficient(Patterns, ...
            'compressReals', true, 'verbose', false);
Patterns = nd.apply(Patterns, "nSource-X_source", @(x) size(x,1) );
Patterns = nd.apply(Patterns, "nTarget-X_target", @(x) size(x,1) );
% Acquire some useful information
uAnimals = unique([Option.animal]);
nSource = zeros(1,nDataset);
nTarget = zeros(1,nDataset);
numDimsUsedForPrediction = cell(1,nDataset);
for a = 1:nDataset
    nSource(a) = size(Patterns(a,1,1,1,1).X_source,1);
    nTarget(a) = size(Patterns(a,1,1,1,1).X_target,1);
    numDimsUsedForPrediction{a} = 1:min(nSource(a), nTarget(a));
end
regions = [const.HPC, const.PFC];

% No idea why sometimes this is fulla nans
if any(ismissing(Option(1).patternNamesFull))
    trans = @(x,y) x + " " + y;
    Option = nd.apply(Option, "patternNamesFull-patternNames,genH_name", trans);
    Option = nd.apply(Option, "patterNamesFull-patternNames,genH_name", trans);
    assert(all(~ismissing(Option(1).patterNamesFull)), "Example number is missing or NaN")
end

% little check
% [[Patterns(:, 1,1,1,1,1,1).animal]',    [Option(:,1,1,1,1,1).animal]']
% [[Patterns(:, 1,1,1,1,1,1).generateH]', [Option(:,1,1,1,1,1).generateH]']
[~,I]  =sortrows([[Option.generateH]; [Option.animal]]');
Option = Option(I);
Patterns = Patterns(I,:,:,:,:,:,:);
% ISSUE: THIS BREAKS FOR MIDPATTERN=1 and ZSCR=1
P = munge.reshapeByLabels_v2(Patterns, 1, ["genH_name","animal"]);
O = munge.reshapeByLabels_v2(Option, 1,   ["genH_name","animal"]);
clear Out

%% Calculate the pattern 
T = query.getPatternTable(Patterns, Option);
!pushover-cli "Finished PreambFigs"
disp("Took " + toc + " seconds to load patterns")
