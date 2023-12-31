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
Option.tableAppend  = "_cohpow";
Option.analysis.cca = true;
Option.midpattern   = true;
Option.saveRaw = true;
Option.quantileToMakeWindows = 0.75;
Option.positiveDerivativeCheck = [0 0];
tableCheck = true; % Set to true if you want to check RunsSummary table
X = datetime('now') - hours(3); % Set to the minimum timestamp you want to check RunsSummary table for
set(0, 'DefaultFigureVisible', 'off');
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
    "fromCP  fromRipTimes"... COHERENCE controlled by power
    ..."fromWpli  fromRipTimes", ...    WPLI
];
zvals = [true, false];


%% Print what we're doing
disp(" ----------- RunAll ----------------------")
disp("Running " + numel(animal_list) + " animals");
disp("with "    + numel(h_methods)   + " methods");
disp(" and "    + numel(zvals)       + " zscore options");
disp("and Option struct ")
disp(rmfield(Option, {'animal', 'generateH', 'genH_name'}));
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
[cntAn, cntH, cntZ]         = deal(0);
for zsc = progress(zvals,'Title','zscore');  cntZ = cntZ + 1;
for genH_= progress(h_methods,'Title','genH method'); cntH = cntH + 1;
for iAnimal = progress(1:numel(animal_list),'Title','Animal'); cntAn = cntAn + 1;
        Option.preProcess_zscore = zsc;
        Option.animal = animal_list(iAnimal);
        Option.generateH = genH_;
        disp("-------------------------------");
        disp("GenerateH = " + Option.generateH);
        disp("-------------------------------");
        if tableCheck
            % Check if the combination of animal and generateH exists in RunsSummary and its timestamp is not older than X
            dates = NaT(numel(RunsSummary.timestamp),1);
            for i = 1:numel(RunsSummary.timestamp)
                dates(i) = datetime(RunsSummary.timestamp(i));
            end
            mask = RunsSummary.animal == Option.animal & RunsSummary.generateH == Option.generateH & RunsSummary.midpattern == Option.midpattern & RunsSummary.preProcess_zscore == Option.preProcess_zscore & dates > X;
            if any(mask)
                disp("")
                disp("Skipping as entry found in RunsSummary with recent timestamp");
                pause(0.25)
                continue
            else
                pause(0.25)
                disp("No entry found in RunsSummary with recent timestamp... running...");
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
        table.combineAndUpdateTables("DetailedRunsSummary_*", "DetailedRunsSummary");
        disp("finished " + Option.animal + " " + Option.generateH + " " + Option.preProcess_zscore + "<---" + string(datetime('now')));
        close all
end % genH
end % animal
end % zscore
