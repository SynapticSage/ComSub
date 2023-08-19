function [combinedPatternsTable, singlePredTable] = allAnimPred(Patterns, Option)
    % Get the size of the original array                                                                                                                                                                                                                                        
    const = option.constants();
    regions = [const.HPC, const.PFC];

    sz = size(Patterns);
    % Calculate the product of the dimensions to collapse                                                                                                                                                                                                                       
    collapsed_size = prod(sz(1:2));                                                                                                                                                                                                                                            
    % Reshape the array
    ReshapedPatterns = squeeze(reshape(Patterns, [collapsed_size, sz(end-2:end)]));
    % Preallocate growing vectors
     % Same code as before until line:
    % combinedPatternsTable   = table();
    combinedPatternsTable   = cell(1, numel(ReshapedPatterns));
    % nPatterns = size(ReshapedPatterns, 1);
    % nRegions = length(regions);
    % nPatternAndControl = Option(1).nPatternAndControl;
    % Get the maximum possible nSource
    % max_nSource = max(arrayfun(@(x) size(x.X_source', 2), ReshapedPatterns(:)));
    % Preallocation using the maximum nSource
    singlePredTable = cell(1, numel(ReshapedPatterns));
    inds = 1:numel(ReshapedPatterns);
    inds = inds(:);

    P = ProgressBar(numel(inds));

    parfor p = inds'
        % Existing code here
        % Get current nTarget and nSource
        currSource = ReshapedPatterns(p).X_source';
        currTarget = ReshapedPatterns(p).X_target';
        currB = ReshapedPatterns(p).rankRegress.B;
        nSource = repmat(size(currSource, 2), size(currTarget, 2), 1);
        nTarget = repmat(size(currTarget, 2), size(currTarget, 2), 1);
        iTarget = 1:nTarget;
        % iSource = 1:nSource;
        animal = repmat(ReshapedPatterns(p).animal, nTarget(1), 1);
        genH   = repmat(ReshapedPatterns(p).genH_name, nTarget(1), 1);
        direction = repmat(ReshapedPatterns(p).directionality, nTarget(1), 1);
        name = repmat(ReshapedPatterns(p).name, nTarget(1), 1);
        [perf, mu, dev] = plots.calculatePredictionPerformance(currSource, currTarget, currB);
        mu = repmat(mu, nTarget(1), 1);
        dev = repmat(dev, nTarget(1), 1);
        iTarget = iTarget';
        perf = perf';
        newtab = table(animal, genH, name, direction, nSource, nTarget, iTarget, perf, mu, dev);
        combinedPatternsTable{p} = newtab;
        singlePredTable{p} = cell(1, nSource(1));
        for is = 1:nSource(1)
            if ~isfield(ReshapedPatterns(p), 'rankRegress') || ...
                ~isfield(ReshapedPatterns(p).rankRegress, 'singlesource_B')
                continue
            end
            curr_singleB = ReshapedPatterns(p).rankRegress.singlesource_B;
            if isempty(curr_singleB)
                continue
            end
            curr_singleB = curr_singleB{is};
            curr_singlesource = currSource(1:size(currSource,1),is);
            [perf, mu, dev] = plots.calculatePredictionPerformance(curr_singlesource, currTarget, curr_singleB);
            mu = repmat(mu, nTarget(1), 1);
            dev = repmat(dev, nTarget(1), 1);
            iSource = repmat(is, nTarget(1), 1);
            perf = perf';
            newtab = table(animal, genH, direction, nSource, nTarget, iSource, iTarget, perf, mu, dev);
            singlePredTable{p}{is} = newtab;
        end
        updateParallel([], pwd);
        % P.step([], [], []);
    end
    P.release();
    combinedPatternsTable = vertcat(combinedPatternsTable{:});
    for i = 1:numel(singlePredTable)
        singlePredTable{i} = vertcat(singlePredTable{i}{:});
    end
    combinedPatternsTable.Var8 = combinedPatternsTable.perf;
    combinedPatternsTable.Var7 = combinedPatternsTable.iTarget;
    singlePredTable       = vertcat(singlePredTable{:});
    if ~isempty(singlePredTable)
        % Just to not break legacy downstream code
        singlePredTable.Var9 = singlePredTable.perf;
        singlePredTable.Var8 = singlePredTable.iTarget;
    end
    if ~exist(figuredefine("tables"), 'dir'), mkdir(figuredefine("tables")); end
    assert(~isempty(combinedPatternsTable), 'No data to save');
    writetable(combinedPatternsTable, fullfile(figuredefine("tables"), 'fig2_prediction.csv'));
    if ~isempty(singlePredTable)
        writetable(singlePredTable, fullfile(figuredefine("tables"), 'fig2_singlePrediction.csv'));
    end
    combinedPatternsTable = readtable(fullfile(figuredefine("tables"), 'fig2_prediction.csv'));
    singlePredTable = readtable(fullfile(figuredefine("tables"), 'fig2_singlePrediction.csv'));
