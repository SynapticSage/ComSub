clear
commsubspaceToPath
% RunAll.m
% Runs TheScript.m for all animals and all methods
%
% Assumptions:
%   - You are in the parent ComSub directory
%   - You have added all the subdirectories to the path
%           e.g. addpath(genpath(parent_directory))

%% Script parameters
Option = option.defaults();
Option.tableAppend = "_coh";

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

%  Load previous progress
progress_file = fullfile(hashdefine(), "last_run_coh.mat");
if exist(progress_file, "file")
    load(progress_file, "last_run");
else
    last_run = [];
end

%% Print what we're doing
disp(" ----------- RunAll ----------------------")
disp("Running " + numel(animal_list) + " animals");
disp("with " + numel(h_methods) + " methods");
disp("and Option struct ")
disp(rmfield(Option, {'animal', 'generateH'}));
disp("and analysis struct ")
disp(Option.analysis);
disp("Last run was " + last_run);
disp(" -----------------------------------------")
disp("Press any key to continue");
pause

%% Run
pickup_where_left_off = true;
for genH_= progress(h_methods,'Title','genH method')
    for iAnimal = progress(1:numel(animal_list),'Title','Animal')
        Option.animal = animal_list(iAnimal);
        Option.generateH = genH_;
        disp("Running " + Option.animal + " " + Option.generateH);
        if isequal([Option.animal, Option.generateH], last_run)
           pickup_where_left_off = true;
        else
           disp("Skipping previous run");
        end
        if pickup_where_left_off == false
            continue
        end
        %try
        TheScript
        %catch MatlabException
            %warning('Failed to run %s %s', Option.animal, Option.generateH);
        %end
        last_run = [Option.animal, Option.generateH];
        save(progress_file, "last_run");
    end 
end