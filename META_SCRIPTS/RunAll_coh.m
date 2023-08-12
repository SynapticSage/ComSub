clear
commsubspaceToPath
% RunAll.m
% Runs TheScript.m for all animals and all methods
%
% Assumptions:
%   - You are in the parent ComSub directory
%   - You have added all the subdirectories to the path
%           e.g. addpath(genpath(parent_directory))

%% Combine tables if they exist
table.combineAndUpdateTables("RunsSummary_*", "RunsSummary");
load("RunsSummary.mat", "RunsSummary");
% table.combineAndUpdateTables("DetailedRunsSummary_*", "DetailedRunsSummary");

%% Script parameters
Option = option.defaults();
Option.tableAppend  = "_coh";
Option.analysis.cca = true;
Option.midpattern   = true;
% Option.analysis.checks = true;

animal_list = [...
    "JS21",...
    "JS15",...
    "JS14",...
    "JS13",...
    "JS17",...
    "ER1",...
    "ZT2"...
];

h_methods = [  ...
    ..."fromSpectra  fromRipTimes"   ...SPECTRAL POWER
    "fromCoherence  fromRipTimes"... COHERENCE
    ..."fromWpli  fromRipTimes", ...    WPLI
];
zvals = [true, false];


%% Print what we're doing
disp(" ----------- RunAll ----------------------")
disp("Running " + numel(animal_list) + " animals");
disp("with "    + numel(h_methods)   + " methods");
disp(" and "    + numel(zvals)       + " zscore options");
disp("and Option struct ")
disp(rmfield(Option, {'animal', 'generateH'}));
disp("and analysis struct ")
disp(Option.analysis);
disp(" -----------------------------------------")
disp("Press any key to continue");
pause
first = false;

dopar = false;
if dopar
    jobs = [];
end

%% Run
X = datetime('now') - hours(3);
dates = NaT(numel(RunsSummary.timestamp),1);
for i = 1:numel(RunsSummary.timestamp)
    dates(i) = datetime(RunsSummary.timestamp(i));
end
[cntAn, cntH, cntZ]         = deal(0);
tableCheck = false; % Set to true if you want to check RunsSummary table
for zsc = progress(zvals,'Title','zscore');  cntZ = cntZ + 1;
for genH_= progress(h_methods,'Title','genH method'); cntH = cntH + 1;
for iAnimal = progress(1:numel(animal_list),'Title','Animal'); cntAn = cntAn + 1;
        Option.preProcess_zscore = zsc;
        Option.animal = animal_list(iAnimal);
        Option.generateH = genH_;
        if tableCheck
            % Check if the combination of animal and generateH exists in RunsSummary and its timestamp is not older than X
            mask = RunsSummary.animal == Option.animal & ... 
                   RunsSummary.generateH == Option.generateH & ...
                   RunsSummary.preProcess_zscore == Option.preProcess_zscore & ...
                   dates > X;
            if any(mask)
                disp("Skipping as entry found in RunsSummary with recent timestamp");
                continue
            end
        end
        disp(newline + "-------------------------------");
        disp("Running " + Option.animal + " " + Option.generateH);
        disp("-------------------------------");
        if ~first
            disp("Press any key to continue");
            pause
            first = true;
        end
        if dopar
            jobs(cntH, cntAn) = batch(TheScript, ...
            'Workspace', {Option}, 'CurrentFolder', pwd);
        else
            %try
            if ~exist(figuredefine("logs"), 'dir')
                mkdir(figuredefine("logs"));
            end
            diary(figuredefine("logs", replace(strjoin([Option.animal, Option.generateH, Option.preProcess_zscore, Option.midpattern], "_"), " ", "") + ".log"));
            TheScript
            diary off
            %catch MatlabException
                %warning('Failed to run %s %s', Option.animal, Option.generateH);
            %end
        end
        table.combineAndUpdateTables("RunsSummary_*", "RunsSummary");
        disp("finished " + Option.animal + " " + Option.generateH + " " + Option.preProcess_zscore + "<---" + string(datetime('now')));
        close all
end % genH
end % animal
end % zscore
