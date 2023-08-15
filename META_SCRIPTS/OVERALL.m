py = pyenv("Version", "~/miniconda3/bin/python")
disp(pyversion); % Check python version

% HOUSES A FLOW TO RECREATE ANALYSES/FIGS
                                                
% o  /    ,---.          |                        
%   /     |---|,---.,---.|    ,   .,---.,---.,---.
%  /      |   ||   |,---||    |   |`---.|---'`---.
% /  o    `   '`   '`---^`---'`---|`---'`---'`---'
%                             `---'               

RunAll % Calls TheScript using various Option struct options
% RunAll_spectra
% RunAll_coh

% RunAll* calls TheScript over various Option struct regimes, 
% individual animals and network pattern types. These can generate
% plots about those specific runs, and save processed data as hashed
% matfiles, a log of the run, and record the hash and all general
% details about the run into a tabular database.
                                       
addpath(hashdefine())
% o  /    ,---.o                         
%   /     |__. .,---..   .,---.,---.,---.
%  /      |    ||   ||   ||    |---'`---.
% /  o    `    ``---|`---'`    `---'`---'
%               `---'                    
% Documentation of the data generated per animal
% (... more to come ... )

PreambFigs; % Setup for figures

% ------------------------------------------------------s
% Mostly knockoff figures from the inspired paper source 
% (some novel)
% -------------------------------------------------------
SemedoPaperFigures_m

% -------------------------------------------------
% Grammer of graphics summaries
% -------------------------------------------------
GrammPlots

% -------------------------------------------------
% Subspace angle 
% -------------------------------------------------
plots.subspace.angle.Run; % do this   !!! 
% Python plots
% py.importlib.import_module(codedefine("Notebooks","python","create_clustergram.py"));
% py.importlib.import_module(codedefine("Notebooks","python","create_graph.py"));
sub = subdir(char(figuredefine("subspaceAngle","*.mat")));
% Remove any matfiles with archive string
sub = sub(~contains({sub.name}, "archive") ...
        & contains({sub.name}, "coh"));
for i = progress(1:length(sub),'Title', 'Subspace angle plots')
    file = sub(i).name;
    disp("--------------------");
    disp("TRYING " + file);
    system("python " + codedefine("Notebooks","python","subspace_angle.py") + "  --file " + file);
end


% -----------------------------------------------------
% Python  data combine
% -----------------------------------------------------

command = sprintf('python %s', codedefine("Notebooks","python","uvevents_combine.py"));
system(command);

command = sprintf('python %s', codedefine("Notebooks","python","ccatime_combinefiles.py"));
system(command);

command = sprintf('python %s', codedefine("Notebooks","python","regress_combine.py"));
system(command);

% -----------------------------------------------------
% Python  figures
% -----------------------------------------------------

command = sprintf('python %s', codedefine("Notebooks","python","dimPred.py"));
system(command);

command = sprintf('python %s', codedefine("Notebooks","python","prediction.py"));
system(command);

command = sprintf('python %s', codedefine("Notebooks","python","ccatime.py"));
system(command);

command = sprintf('python %s', codedefine("Notebooks","python","uvevents_plot.py"));
system(command);
% ------------------------------------------------------


tic; load RunsSummary.mat; toc
tic; load DetailedRunsSummary.mat; toc
