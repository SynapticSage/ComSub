% Connecting ComSub with Rhythms 🗣️🔗🧠🎼
% ========================================
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
    disp("Option struct does not exist, creating one")
    addpath(genpath(codedefine()));
    Option = option.defaults(); 
else
    disp("Option struct already exists, using that")
    Option = option.setdefaults(Option);
    disp("Option struct is: ")
    disp(Option)
end
hash = store.gethash(Option);
disp("Hash is: " + hash)
if ~exist('RunsSummary','var')
    load RunsSummary
    if ismember(hash, RunsSummary.hash), disp("config exists!")
    else, disp("new configuration")
    end
end
%%%%%%%%%%%%%%%% DISPLAY OUR OPTIONS TO USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(Option.loadifexists, false) && ...
    exist(hash + ".mat", 'file')
    disp("Loading from file: " + store.gethash(Option) + ".mat")
    m = matfile(store.gethash(Option) + ".mat");
    % m = matfile("bef0923.mat", "Writable", true);
    disp("Loaded variables: ")
    Events             = util.matfile.getdefault(m, 'Events', []);
    Spk                = util.matfile.getdefault(m, 'Spk', []);
    Patterns           = util.matfile.getdefault(m, 'Patterns', []);
    Patterns_overall   = util.matfile.getdefault(m, 'Patterns_overall', []);
    Components         = util.matfile.getdefault(m, 'Components', []);
    Components_overall = util.matfile.getdefault(m, 'Components_overall', []);
    Option             = util.matfile.getdefault(m, 'Option', []);
    disp("...done")
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
    Components = nd.initFrom(Patterns, ...
    {'index_source', 'index_target', 'directionality', 'name'});
    Components_overall = nd.initFrom(Patterns_overall, ...
    {'index_source', 'index_target', 'directionality', 'name'});
end

%%%%%%%%%%%%%%%% ANALYSIS SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------")
disp("     Analysis           ")
disp("------------------------")

if Option.analysis.rankRegress % 
    % Rank regression of network pattern windows of spiking activity
    % (Subspaces acquired here)
    % TODO: 
    % 1. fix Option.rankregress => Option.rankRegress                 
    % 2. most rankRegress.B_ are empty                                
    Patterns         = analysis.rankRegress(Patterns, Option, 'verbose', true);        
    Patterns_overall = analysis.rankRegress(Patterns_overall, Option);
end

if Option.analysis.factorAnalysis % 󰍉
    % Factor analysis of network pattern windows of spiking activity
    % (Used to measure instrinsic dimensionality of network activity)
    Patterns = analysis.factorAnalysis(Patterns, Option);
end

if Option.analysis.cca % 🔀
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    % Some metadata for steps below
    figAppend    = strjoin([Option.animal,Option.genH_name,Option.preProcess_zscore], "_");
    structAppend = struct('animal', Option.animal, 'genH', Option.genH_name, 'zscore', Option.preProcess_zscore);
    % ---------------------------------------------------------------------
    % (Subspaces acquired here)
    % Patterns         = analysis.cca(Patterns, Option);
    Patterns_overall = analysis.cca(Patterns_overall, Option);
    % Copy over the components
    Components_overall =  ...
         nd.fieldSet(Components_overall, 'cca', Patterns_overall);
    analysis.cca.savedat(Patterns_overall, Spk, behavior, Option, figAppend);
    % ---------------------------------------------------------------------
    % Event analysis ------------------------------------------------------
    % (append cca commsub levels during events)
    % 藺
    event_anal  = ... 
         analysis.cca.event_analysis(Patterns_overall, Spk, Events, Option, behavior);
    for i = 1:Option.nPatternAndControl+1
        plots.plot_event_values(event_anal(2,i), 'figAppend', figAppend + "_i");
    end
    Components_overall = ... 
         nd.fieldSet(Components_overall, 'event_anal', event_anal);
    % ---------------------------------------------------------------------
    % Create table of results ---------------------------------------------
    % (create a table regarding cca versus efizz and behavior)
    % 藺 ccatime
    efizz = load(Option.animal + "spectralBehavior.mat", "efizz");
    efizz = efizz.efizz;
    table.analyses.ccatime(Patterns_overall, Spk, efizz, Option, behavior,...
                          'behaviorColumns', ...
    {'vel', 'accel', 'lindist', 'rewarded', ...
    'trajbound','inBoundChoiceTimes','outBoundChoiceTimes','rewardTimes'});
    % ---------------------------------------------------------------------
    % Triggered spectrogram u----------------------------------------------
    % (create compute triggered spectrograms for commsubs)
    disp("Running triggered spectrogram - run")
    close all
    % 藺 triggered win
    triggered_spectrogram_run = ...
    analysis.cca.triggered_spectrogram(Patterns_overall, Spk, efizz,...
    'ploton', true, ... 
    'quantile_threshold', 0.85,...
    'figAppend', structAppend, ...
    'based_on', 'mean_and_std', ...
    'runtype', 1);
    if ~exist(figuredefine("data"), 'dir'); mkdir(figuredefine("data")); end
    save(figuredefine("data", "trigspec_" + figAppend), "Option", "triggered_spectrogram_run");
    dcnt=0;
    % 藺 regress
    tic;
    for d = progress([inf], 'Title', 'Regress-faxis'); dcnt=dcnt+1;
        for i = progress(1:size(Patterns_overall,2), 'Title', 'Regress')
            for f = progress(["phi","S1","S2","Cavg","wpli_avg"],'Title', 'Regress-field')
                Patterns_overall(2,i).regress(dcnt).(f) = ... 
                analysis.cca.regressefizz(efizz, Patterns_overall(2,i), f,...
                'faxis', d, "tabPrepend", figAppend, 'ploton', true);
            end
        end
    end
    close all
    disp("Regressed efizz in " + toc + " seconds")
    % ---------------------------------------------------------------------
end

if Option.analysis.timeVarying % 🎶
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    % How much spiking moment to moment is explained by subspace
    % Requirements: Option.analysis.cca
    % TODO: 1 .also return epochwise zscored neural firing matching
    %       2. return timeseries of smoothed firing rate
    % Component matching over time
    rrr = analysis.match_rrr(Patterns_overall, Option, Spk);
    cca = analysis.match_cca(Patterns_overall, Option, Spk);
    % Spectral matches
    Components_overall = ... 
    plots.temporal.correlateSpectral(Components_overall, Events, Option);
    Components_overall = ... thescript
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
    else; wait_state = false; end
    plots.runChecks(Events, Spk, Patterns, Option, ...
                    'parallel', false, 'wait', wait_state, 'visible', 'off');
end

% TODO: (1) plug in JPECC version of rankRegress here
% TODO: (2) function that outputs average response of Pattern struct per neuron
close all

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Option.save % 💾

    disp("Saving results")
    store.savetables(Events, Patterns, Option);

    saveVars = [];
    if exist('Patterns','var')
        saveVars.Patterns = Patterns;
        saveVars.Patterns_overall = Patterns_overall;
    end
    if exist('Components', 'var')
        saveVars.Components         = Components;
        saveVars.Components_overall = Components_overall;
    end
    store.savevars(Option, Events, Spk, saveVars);
    disp("...done")
    !pushover-cli "finished saving"

end
