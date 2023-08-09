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

% -----------------------------
% Grammer of graphics summaries
% -----------------------------
GrammPlots

% -----------------------------
% Subspace angle 
% -----------------------------
plots.subspace.angle.Run; % do this   !!! 
% Python plots
py.importlib.import_module(codedefine("Notebooks","python","create_clustergram.py"));
py.importlib.import_module(codedefine("Notebooks","python","create_graph.py"));

% --------------
% Python  figures
% --------------
script1 = ...
py.importlib.import_module(codedefine("Notebooks","python","dimPred.py"));
script2 = ...
py.importlib.import_module(codedefine("Notebooks","python","prediction.py"));
script3 = ...
py.importlib.import_module(codedefine("Notebooks","python","ccatime.py"));
script4 = ...
py.importlib.import_module(codedefine("Notebooks","python","uvevents_combine.py"));
script5 = ...
py.importlib.import_module(codedefine("Notebooks","python","uvevents_plot.py"));
% ----------------------------------------------


tic; load RunsSummary.mat; toc
tic; load DetailedRunsSummary.mat; toc
