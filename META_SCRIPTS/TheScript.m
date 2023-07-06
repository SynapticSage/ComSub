% TODO:
% 1. have way of passing generate_H cell, and then just load all of those
%    patterns simulateously in 1 go instead of as separate runs
% 
% ---- PATH -----
% Matlab uses startup.m to run startup code...
% Put this in your startup.m so that the code for this is in path:
%
% addpath(genpath('/Volumes/MATLAB-Drive/')) % or wherever your CODE files are
% located
% addpath(genpath('~/Data/Raw/')) % or wherever your DATA files are located
% ===================================================
% OPTION STRUCT encoding properties of the script run
% ===================================================
% see +option.default() to set default options
if ~exist('Option','var')
    Option = option.defaults(); 
else
    Option = option.setdefaults(Option);
end
%%%%%%%%%%%%%%%% DISPLAY OUR OPTIONS TO USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(Option.loadifexists, false) && exist(store.gethash(Option) + ".mat", 'file')
    disp("Loading from file: " + store.gethash(Option) + ".mat")
    % m = matfile(store.gethash(Option) + ".mat");
    m = matfile("bef0923.mat", "Writable", true);
    disp("Loaded variables: ")
    Events             = util.matfile.getdefault(m, 'Events', []);
    Spk                = util.matfile.getdefault(m, 'Spk', []);
    Patterns           = util.matfile.getdefault(m, 'Patterns', []);
    Patterns_overall   = util.matfile.getdefault(m, 'Patterns_overall', []);
    Components         = util.matfile.getdefault(m, 'Components', []);
    Components_overall = util.matfile.getdefault(m, 'Components_overall', []);
    Option             = util.matfile.getdefault(m, 'Option', []);
    disp("...done")
    clear m
else
    %%%%%%%%%%%%%%%% OBTAIN EVENT MATRICES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("------------------------")
    disp("    Obtaining events    ")
    disp("------------------------")
    Events = events.ThetaDeltaRipple(Option);
    % Documentation
    % Events is a struct with fields:
    % - .times : array of times of events
    % - .H     : Event Matrix,    T x 3, and each column are theta, delta, ripple
    % - .Hvals : Event Matrix,    T x 3, values without nans
    % - .Hnanlocs : Event Matrix, T x 3, logicals of nans

    %%%%%%%%%%%%%%%% CUT WINDOWS WITH EVENT MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("------------------------")
    disp("    Cutting windows     ")
    disp("------------------------")
    Events = windows.ThetaDeltaRipple(Events, Option);
    % -  cutoffs:       nPatterns x 1 vector of cutoffs
    % TODO: modify to be able to include overall pattern and track patterns
    % PRIORITY; overall: medium, track: very low, overall can be included in
    % cellOfWindows, whereas, track can be included as a separate output

    %%%%%%%%%%%%%%%% ACQUIRE SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting spikes
    disp("------------------------")
    disp("    Getting spikes      ")
    disp("------------------------")
    Spk = spikes.getSpikeTrain(Option.animal, Option.spikeBinSize, ...
                             Option.samplingRate);

    % filter the neurons whose firing rate is lower than specified threshold
    if Option.preProcess_FilterLowFR 
        disp("------------------------")
        disp("Filtering low FR neurons")
        disp("------------------------")
        Spk = trialSpikes.filterFR(Spk, 0.1);
        disp("Mean FR: " + sort(Spk.avgFR))
    end

    if Option.preProcess_gaussianFilter
        % Gaussian filter the spikeCountMatrix/spikeRateMatrix
        gauss = gausswin(Option.preProcess_gaussianFilter);
        for i = progress(1:size(Spk.spikeRateMatrix, 1), 'Title', 'Gaussian filtering')
            Spk.spikeRateMatrix(i, :)  = conv(Spk.spikeRateMatrix(i, :), gauss, 'same');
            Spk.spikeCountMatrix(i, :) = conv(Spk.spikeCountMatrix(i, :), gauss, 'same');
        end
    end

    if Option.preProcess_zscore
        % Z-score the spikeCountMatrix/spikeRateMatrix
        disp(" Z-scoring ")
        if ~isfield(Spk, 'muFR')
            Spk.muFR  = mean(Spk.spikeRateMatrix, 2);
            Spk.stdFR = std(Spk.spikeRateMatrix, 0, 2);
        end
        Spk.spikeRateMatrix  = zscore(Spk.spikeRateMatrix,  0, 2);
        Spk.spikeCountMatrix = zscore(Spk.spikeCountMatrix, 0, 2);
        Spk.avgFR = mean(Spk.spikeRateMatrix, 2);
    end
    prewindow_copy = Spk;

    % %%%%%%%%%%%%%% ACQUIRE TRIALS FROM WINDOWS + SPIKES %%%%%%%%%%%%%%%%%%%
    % RYAN bug here .. timeBinStartEnd instead of timeBinMidPoints
    disp("------------------------")
    disp("   Windowing spikes     ")
    disp("------------------------")
    Spk = trialSpikes.generate(Spk, Events, Option);

    %%%%%%%%%%%%%%%%% SETUP RAW DATA STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Structure for separated data
    %%%%%%%%%%%%%%%% SEPRATE BRAIN AREA DATASULT STRUCTURES %%%%%%%%%%%%%%%%%%
    % Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and
    % neurons that in HPC
    %% Separate firing pattern into source and target
    [Spk.nSource,~,~] = size(Spk.hpc.X{1});
    [Spk.nTarget,~,~] = size(Spk.pfc.X{1});
    Spk.celllookup = cellInfo.getCellIdentities(Option.animal, Spk.cell_index,...
                                                Spk.areaPerNeuron);
    system("pushover-cli 'Finished munging data for analysis'");

    %%%%%%%%%%%%%%%% SETUP PARTITIONS AND RESULT STRUCTURES %%%%%%%%%%%%%%%%%%
    disp("------------------------")
    disp(" Subsampling partitions ")
    disp("------------------------")
    [Patterns, Patterns_overall] = trialSpikes.partitionAndInitialize(Spk, Option);
end

%%%%%%%%%%%%%%%% ANALYSIS SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------")
disp("     Analysis           ")
disp("------------------------")

if Option.analysis.rankRegress
    % Rank regression of network pattern windows of spiking activity
    % (Subspaces acquired here)
    % TODO: 
    % 1. fix Option.rankregress => Option.rankRegress                 
    % 2. most rankRegress.B_ are empty                                
    Patterns         = analysis.rankRegress(Patterns, Option);        
    Patterns_overall = analysis.rankRegress(Patterns_overall, Option);
end

if Option.analysis.factorAnalysis
    % Factor analysis of network pattern windows of spiking activity
    % (Used to measure instrinsic dimensionality of network activity)
    Patterns = analysis.factorAnalysis(Patterns, Option);
end

if Option.analysis.cca
    % ISSUE: warnings emitted regarding full rank -- linear independence
    % violation can be subtle problem or not a problem at all
    % remedies (1) regularize (2) remove linearly dependent columns (PCA)
    Patterns         = analysis.cca(Patterns, Option);
    Patterns_overall = analysis.cca(Patterns_overall, Option);
    for i = 1:numel(Patterns_overall)
        Components_overall(i).cca = Patterns_overall(i).cca;
    end
    % TODO : section that knocks off kim 2022 after these measurements
    Components.cca   = analysis.cca.event_analysis(Patterns_overall, Spk, Events, Option);
end

if Option.analysis.timeVarying
    % How much spiking moment to moment is explained by subspace
    % ISSUE: hits a bug on line 4
    % TODO: 1 .also return epochwise zscored neural firing matching
    %       2. return timeseries of smoothed firing rate
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    % Component matching over time
    rrr = analysis.match_rrr(Patterns_overall, Option, Spk);
    cca = analysis.match_cca(Patterns_overall, Option, Spk);
    % Spectral matches
    Components_overall = ... 
    plots.temporal.correlateSpectral(Components_overall, Events, Option);
    Components_overall = ... 
    plots.temporal.correlateSpectral(Components_overall, Events, Option, 'componentMethod', 'cca');
    % Behavior matches
    Components         = plots.temporal.correlateBehavior(Components, Events, Option);
end

if Option.analysis.checks
    % Plots regarding the raw and processed data (and sometimes
    % relation to processed Patterns struct)
    % TODO: Think about splitting this into checks involving
    %        versus not involving the Patterns struct
    if strcmp(Option.animal, "JS21"); wait_state = true;
    else; wait_state = false; end;
    plots.runChecks(Events, Spk, Patterns, Option, ...
                    'parallel', true, 'wait', wait_state);
end

% TODO: (1) plug in JPECC version of rankRegress here
% TODO: (2) function that outputs average response of Pattern struct per neuron

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Option.save

    disp("Saving results")
    store.savetables(Events, Patterns, Option);

    saveVars = {};
    if exist('Patterns','var')
        saveVars = [saveVars, {'Patterns', Patterns, ...
                               'Patterns_overall', Patterns_overall}];
    end
    if exist('Components', 'var')
        saveVars = [saveVars, {'Components', Components, ...
                               'Components_overall', Components_overall}];
    end
    store.savevars(Option, Event, Spk, saveVars{:});
    
end
