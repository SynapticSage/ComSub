%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Mode %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
Option.waysOfPartitions = 2; %whether the script is bidirectional (split neurons in two ways)
run_selected_genH = false;

%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%% Paths %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('pathSet','var') || pathSet == 0 % only set path if unset
    addpath(genpath(pwd)) % all folders in utils added, including semedo code
    
    % Whose computer are we running?
    if ispc
        paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
    elseif ismac
        paths(1) = "/Volumes/sharespace-commsub/data";
        paths(2) = "~/Data/commsubspace";
    else
	paths = [];
    end
    % Set data paths for that computer
    arrayfun(@(path) addpath(genpath(path)), paths);
    % ---------- Add paths for specific users ------------------------------
    % addpath(genpath(datadef(user))); % this line throws an error for pc
    
    addpath(...
        fullfile(codedefine,"Shared"));
    addpath(genpath(fullfile(codedefine, "Shared", "utils")));
    pathSet = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%
%% Script parameters %%
%%%%%%%%%%%%%%%%%%%%%%%
Default = struct();
Default.animal = "ZT2";
% determined by spectral power with pattern determined by global ripple
% Default.generateH = "fromSpectra "+" fromRipTimes";
% Default.generateH = "fromWpli " + " fromRipTimes";
% Default.generateH = "fromCoherence "+" fromRipTimes";
Default.samplingRate  = [] ;              % For spikes.getSpikeTrain, nan if not given
Default.spikeBinSize  = 0.15;               % 100 milliseconds
Default.timesPerTrial = 10;                % 10 times per trial
Default.winSize       = {[-0.10, 0.20]};     % size of the window
Default.sourceArea    = "CA1";             % only when there are
Default.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
Default.singleControl = false;                 % whether to use just one control column
Default.numPartition = 200;                    % ways to split source and target
Default.binsToMatchFR = 20;
Default.lowerControl = true;
Default.preProcess_FilterLowFR = true;
Default.preProcess_matchingDiscreteFR = true;
Default.oldControlBehavior = false;
Default.coherence_pattern = "wpli";

% %%%%%%%%%%%%%%%%%%%%%%%
% %% Set offical option %
% %%%%%%%%%%%%%%%%%%%%%%%
% % Write defaults to option struct if not provided
% if ~exist('Option','var')
%     Option = Default;
% else
%     for field = string(fieldnames(Default))'
%         if ~isfield(Option, field)
%             Option.(field) = Default.(field);
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Aliases and Shortcuts %
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% animal = Option.animal;
% animalToPath(animal);
% disp("Running with Option struct => ")
% disp(Option);
% %% Shortcut/alias variables to improve readability
% THETA  = 1;
% DELTA  = 2;
% RIPPLE = 3;
% if Option.sourceArea == "CA1"
%     HPC = 1;
%     PFC = 2;
% else
%     PFC = 1;
%     HPC = 2;
% end
% patternNames = ["theta","delta","ripple"];
% 
% frequenciesPerPattern = [6 10; 0.5 4; 150 200];
% [nPatterns,~] = size(frequenciesPerPattern);
% 
% if Option.singleControl == true
%     Option.nPatternAndControl = nPatterns+1;
% else
%     Option.nPatternAndControl = nPatterns*2;
% end
% patternNames = ["theta","delta","ripple",...
%                 "theta-control","delta-control","ripple-control"];
% if Option.nPatternAndControl == nPatterns+1
%     patternNames = ["theta","delta","ripple","control"];
% end
% 
% animal_list = ["JS21", "ZT2","JS15","JS14","ER1"];
% for animal = animal_list
%     animalToPath(animal);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Compute our pattern-specific spiking data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% patternNames = patternNames(1:3);
% Option.quantileResolution = 3;
% quantile_edges = linspace(0.01,0.99,Option.quantileResolution+1);
% quantile_edges = [quantile_edges(1:end-1)', quantile_edges(2:end)'];
% clear overall;
% for animal = progress(animal_list,'Title','Animals')
%     
%     animalToPath(animal);
%     load(animal + "avgeeg.mat");
%     load(animal + "spectralBehavior.mat");
%     load(animal + "globalripple01.mat");
%     [~, sleepSessions] = getRunningSessions(animal);
% 
%     overall.(animal).cellOfWindows = {};
%     overall.(animal).cutoffs = {};
%     for hilbertPower = progress(1:Option.quantileResolution,'Title','hilbert quantiles')
%         for coherence = 1:Option.quantileResolution
% 
%             p = 0;
%             for pattern = patternNames(1:3); p = p+1;
% 
%                 % ----
%                 % Hilb
%                 % ----
%                 [Hstruct.(pattern).hilb.H, Hstruct.(pattern).hilb.Hvals, ...
%                  Hstruct.(pattern).hilb.Hnanlocs, Hstruct.(pattern).hilb.Htimes] =...
%                     eventMatrix.generateFromFilteredEEG(avgeeg, Option.sourceArea, "patterns", pattern, "downsample", 10, "sessionsToExclude", sleepSessions);
% 
%                 %---
%                 %Coh
%                 %---
%                 if Option.coherence_pattern == "coherence"
%                     spectrogram = efizz.C;
%                 elseif Option.coherence_pattern == "wpli"
%                     spectrogram = efizz.wpli;
%                 else
%                     error("Core method for deriving the event matrix is not recognized");
%                 end
%                 Htimes = efizz.t;
%                 frequencyAxis = efizz.f;
%                 [Hstruct.(pattern).coh.H, Hstruct.(pattern).coh.Hvals, ...
%                 Hstruct.(pattern).coh.Hnanlocs, Hstruct.(pattern).coh.Htimes] = ...
%                     eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
%                     frequenciesPerPattern(p,:));
% 
%                 if pattern == "ripple"
% 
%                     %% Modify ripple pattern? Lower the threshold as with window siz
%                     rip_kws = {'amplitude_at_riptime', true};
% 
%                     [Hstruct.(pattern).coh.Htimes, Hstruct.(pattern).coh.H, ...
%                      Hstruct.(pattern).coh.Hnanlocs, Hstruct.(pattern).coh.Hvals, ...
%                      Hstruct.(pattern).coh.minRippleThreshold, original] = ...
%                         eventMatrix.generateFromRipples(globalripple, ...
%                         rip_kws{:}, ...
%                         'globalrippleWindowUnits', 'amp', ...
%                         'rippleBandTime', Hstruct.(pattern).coh.Htimes,...
%                         'rippleBand', Hstruct.(pattern).coh.Hvals); 
% 
%                     [~, Hstruct.(pattern).hilb.H, ...
%                      Hstruct.(pattern).hilb.Hnanlocs, Hstruct.(pattern).hilb.Hvals, ...
%                      Hstruct.(pattern).hilb.minRippleThreshold, original] = ...
%                         eventMatrix.generateFromRipples(globalripple, ...
%                         rip_kws{:}, ...
%                         'rippleBandTime', Hstruct.(pattern).hilb.Htimes,...
%                         'rippleBand', Hstruct.(pattern).hilb.Hvals,... %RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
%                         'globalrippleWindowUnits', 'std');
% 
%                     % Apply special ripple windowing settings?
%                     Hstruct.(pattern).hilb.arrayForQuantile = Hstruct.(pattern).hilb.H; % Can only use actual ripple times
%                     Hstruct.(pattern).coh.arrayForQuantile  = Hstruct.(pattern).coh.Hvals;
%                 end
%             end
% 
%             %% -----------------------
%             %% Generate joint windows!
%             %% -----------------------
%             % Encode settings for making windows
%             ipat = 0;
%             for pattern = patternNames(1:3)
%                 ipat = ipat + 1;
%                 Hstruct.(pattern).hilb.quantileBin = quantile_edges(hilbertPower, :);
%                 Hstruct.(pattern).coh.quantileBin  = quantile_edges(coherence, :);
%             end
%             % Generate windows
%             [cellOfWindows, cutoffs] = windows.jointMake(Hstruct, Option.winSize, 'patternNames', patternNames);
%             nWindows = cellfun(@(x) size(x,1), cellOfWindows);
%             disp(join(patternNames + " " + nWindows, newline));
% 
%             overall.(animal).cellOfWindows(hilbertPower, coherence, :) = cellOfWindows;
%             overall.(animal).cutoffs{hilbertPower, coherence}          = cutoffs;
%             overall.(animal).nWindows(hilbertPower, coherence,:)       = nWindows(:);
%         end
%     end
%     overall.(animal).cutoffMat = cell2mat(cellfun(@(x) permute(x, [3 4 1 2]), overall.(animal).cutoffs, 'UniformOutput', false));
%     overall.(animal).cutoffMat = permute(overall.(animal).cutoffMat, [3 1 2 4]);
% 
% end
%%
fig('Windows across cond');clf;
t=tiledlayout(3,5, 'TileIndexing', 'columnmajor');
for animal = animal_list
    for i = 1:3; 
        nexttile; 
        imagesc(log10(overall.(animal).nWindows(:,:,i)+10)); 
        xlabel("Power level");
        ylabel("Coherence Level");
        %xticks(1:10);yticks(1:10);
        %XL = num2cell([overall.cutoffMat(:, 1, 1, 1), overall.cutoffMat(:, 1, 1, 2)],2);
        %XL = cellfun(@(x) char("[" + join(string(x),",") + "]"), XL,'UniformOutput',false);
        %YL = num2cell([squeeze(overall.cutoffMat(1, :, 1, 1)), squeeze(overall.cutoffMat(1, :, 1, 2))],2);
        %YL = cellfun(@(x) char("[" + join(string(x),",") + "]"), YL,'UniformOutput',false);
        %xticklabels();yticklabels(1:10);
        axis xy
        title(animal + newline+ patternNames(i)); colorbar; 
    end; 
    sgtitle('Number of windows');
    % (EVENTUALLY, WE WILL TOSS ANY COMBINATIONS WITH WINDOWS LESS THAN RANK(B))
end


%%
Raw = struct();


%% Within each animal balance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% SPIKE SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for animal = animal_list(1:end)
    %% Getting spikes
    [r.timeBinStartEnd, r.timeBinMidPoints, ~, r.spikeCountMatrix, r.spikeRateMatrix, ...
        r.areaPerNeuron, r.cell_index, r.sessionTypePerBin] = spikes.getSpikeTrain(Option.animal, ...
        Option.spikeBinSize,  Option.samplingRate);

    if Option.preProcess_FilterLowFR % filter the neurons whose firing rate is lower than specified threshold
        [r.spikeCountMatrix, r.spikeRateMatrix, r.avgFR, r.areaPerNeuron, r.cell_index]...
            = trialSpikes.filterFR(r.spikeCountMatrix, r.spikeRateMatrix, 0.1, ...
                                   r.timeBinStartEnd, r.areaPerNeuron, r.cell_index);
    end
    Raw.(animal) = r;
end

%% Balance windows
for animal = animal_list
    
    r = Raw.(animal);
    % Balance cecll of windows
    cellCounts = arrayfun(@(x) sum(r.areaPerNeuron==x), unique(r.areaPerNeuron));
    mCellCount = min(cellCounts);
    overall.(animal).minCellCountThresh = mCellCount;

    % Constrain cellOfWindows : TODO more sophisticated version of this
    cellOfWindows = overall.(animal).cellOfWindows;
    nWindows = cellfun(@(x)size(x,1), overall.(animal).cellOfWindows);   %TODO change this to more sophisticated method
    inds = find(nWindows < mCellCount);
    dims = cell(1, ndims(cellOfWindows));
    [dims{:}] = ind2sub(size(cellOfWindows), inds);
    for i = 1:numel(dims{1})
        D = cellfun(@(x) x(i), dims, 'UniformOutput', false);
        cellOfWindows{D{:}} = [];
        nWindows(D{:}) = 0;
    end
    overall.(animal).threshCellOfWindows = cellOfWindows;
    overall.(animal).threshNWindows      = nWindows; 
end

%% OBtain trial windows of neural activity
for animal = animal_list
    r = Raw.(animal);
    %% Cut trials from windows
    [r.spikeSampleMatrix, r.spikeSampleTensor, r.trialTimes] = trialSpikes.generate(...
        r.spikeCountMatrix,...
        r.timeBinStartEnd, overall.(animal).cellOfWindows, ...
        Option.timesPerTrial, Option.nPatternAndControl);
    Raw.(animal) = r;
    
end

%% 
for animal = animal_list
    %p = Patterns.(animal);
    r = Raw.(animal);
    %% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
    [pfc.FR, hpc.FR] = trialSpikes.separateFiringRate(r.avgFR, r.areaPerNeuron);
    pfc.X = trialSpikes.separateSpikes(r.spikeSampleMatrix, r.areaPerNeuron, "PFC");
    hpc.X = trialSpikes.separateSpikes(r.spikeSampleMatrix, r.areaPerNeuron, "CA1");

    %% Separate firing pattern into source and target
    [pfc.nNeurons,~] = size(pfc.X{1});
    [hpc.nNeurons,~] = size(hpc.X{1});
    r.celllookup = cellInfo.getCellIdentities(animal, r.cell_index, r.areaPerNeuron);
    r.pfc = pfc;
    r.hpc = hpc;
    % Add window information to r
    r.windowInfo = overall.(animal);

    %% Partition the cells %%
    p = trialSpikes.partition(r, Option);

    Raw.(animal)      = r;
    Patterns.(animal) = p;
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% RANK REGRESS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for animal = progress(animal_list,'Title','Animal rank regress')

    pat = Patterns.(animal)(:,:,:,:,:); % Focus on cross area interactions
    nTarget = numel(pat(1).index_target);
    nSource = numel(pat(1).index_source);

    if Option.waysOfPartitions ~= 2
        nTarget = size(pat(1).X_target,1);
        nSource = min(size(pat(1,1,1,1,1,1,1).X_source,1),...
            size(pat(1,2,1,1,1,1).X_source,1));
    end
    numDimsUsedForPrediction = 1:min(nTarget,nSource);
    B_singleprediction = cell(1,nSource);
    dim_singleprediction = cell(1,nSource);
    cvNumFolds = 10;

    for p = progress(9626:numel(pat),'Title',char(animal + "rank regress"))


        % Number of cross validation folds.
        cvOptions = statset('crossval');
        regressMethod = @ReducedRankRegress;
        cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
            (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
            numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
            false, 'Scale', false);

        if isempty(pat(p).X_source)
            continue
        end

        % when the partition is three-ways, j==1 means same target/source
        % pair and j==2 means diff target/source pair
        curr_source = (pat(p).X_source)';
        curr_target = (pat(p).X_target)';
        [   pat(p).rankRegress.cvl, ...
            pat(p).rankRegress.cvLoss,...
            pat(p).rankRegress.optDimReducedRankRegress,...
            pat(p).rankRegress.B,...
            pat(p).rankRegress.B_,...
            pat(p).rankRegress.V] ...
            = rankRegressRoutine(cvFun, cvNumFolds, ...
            cvOptions, curr_target, curr_source, ...
            numDimsUsedForPrediction);

        pat(p).X_source = []; % delete to open up space
        pat(p).X_target = [];

        % Single neuron prediction
        %for k = 1:nSource
        %    %curr_singlesource = curr_source(:,j); % RY BUG ALERT!
        %    curr_singlesource = curr_source(:,k); % Fix
        %    if clean.zeroFiring(curr_singlesource)
        %        continue;
        %    end
        %    [~,~, ...
        %        dim_singleprediction{k}, ...
        %        B_singleprediction{k},~,~] = ...
        %        rankRegressRoutine(cvFun, cvNumFolds, ...
        %        numDimsUsedForPrediction);
        %end
        %pat(p).rankRegress.singlesource_B = B_singleprediction;
        %pat(p).rankRegress.singlesource_optDim = ...
        %    dim_singleprediction;
        pat(p).rankRegress.B_rrr = getReducedB_(pat(p).rankRegress.B,...
            pat(p).rankRegress.V, nSource, nTarget,...
            pat(p).rankRegress.optDimReducedRankRegress);
        Patterns.(animal) = pat;
    end

    % Assign performance while removing dims
    %for p = 1:Option.numPartition
    %    try
    %        [pat(p).rankRegress.removedPerformance, ~,~]=...
    %            plots.getUncorrelatedPerformance( pat(p).rankRegress.B_,...
    %            pat(p).X_source, pat(p).X_target, pat(p).rankRegress.optDimReducedRankRegress,...
    %            numDimsUsedForPrediction, pat(p).rankRegress.cvLoss);
    %    catch
    %        keyboard
    %    end
    %    [pat(p).perfAsRemovingEach] = plots.sequentialRemovePredDims(pat(p).X_source, ...
    %        pat(p).X_target,pat(p).rankRegress.B_, pat(p).rankRegress.optDimReducedRankRegress, ...
    %        pat(p).rankRegress.cvLoss,numDimsUsedForPrediction);
    %end
end

%%
% -----|
% ------|
% Plots | 
% -----|
% -----|
fig("Dimension")
for animal = animal_list
    pat = Patterns.(animal);
    indices = nd.indicesMatrixForm(pat);
    D = zeros(size(pat));
    for ind = indices'
        ind = num2cell(ind);
        rr = pat(ind{:}).rankRegress;
        if ~isempty(rr)
            D(ind{:}) = rr.optDimReducedRankRegress;
        end
    end

    reshape_for_plotting = @(x) permute(x, [5 2 3 4 1]); % permute to netpatt direct power coh sample
    D = reshape_for_plotting(D);


    fig("ANimal Dimension")
    for name = 1:numel(patternNames)
    for direct = [1]
        nexttile;
        d = squeeze(D(name, direct, :, :, :));
        x = 1:size(d,1);
        y = 1:size(d,2);
        [X, Y] = ndgrid(x,y);
        Xscat = repmat(X, [1 1 size(d,3)]);
        Yscat = repmat(Y, [1 1 size(d,3)]);
        mu_d = mean(d,3);
        %scatter3(Xscat(:),Yscat(:),d(:))
        %hold on;
        %surf(X, Y, mu_d)
        imagesc(x,y,mu_d)
        axis xy
        xlabel("power")
        ylabel("coherence")
    end
    end


end
