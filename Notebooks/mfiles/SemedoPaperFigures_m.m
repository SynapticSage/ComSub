PreambFigs;

% TODO:
% 1. Animal-wise version of this 
% 2. Add sig to figs
dump = matfile(fullfile(figuredefine("SPF")), "Writable", true);

%%  
% The different animals loaded will actually collapse into the partition
% 
% dimenion, segregated by the genH methods
% 
% *(nMethods * nPartiton * nDirection * nPatterns)*


% Figure 2A: Cofiring in source and target
% 
% Figure 2B: Prediction Performance
% 
% Figure 4A/B: Results from rank regress
% 
% Figure 4C: How optimal number of predicitive dimensions compare between hpc 
% and pfc
% 
% Figure 6: removing predicitive dimensions from both sources
% 
% Figure 7: Predictive performance from dominant dimensions
%% Figure 2
% 
% A: How pairs of neuron in source/target co-fire
% This is method insensitive. Just lump over animals

%% 2A calculate co-firing across all animals
Fig2 = struct();
Fig2.a.all = plots.cf.cofiring(Patterns, Option,...
'appendAttributes', {'animal', 'generateH'});
[spec, opt] = munge.getH(Patterns, Option, "spec");
Fig2.a.spec = plots.cf.cofiring(spec, opt,...
'appendAttributes', {'animal', 'generateH'});
[coh, opt] = munge.getH(Patterns, Option, "coh");
Fig2.a.coh = plots.cf.cofiring(coh, opt);
[wpli, opt] = munge.getH(Patterns, Option, "wpli");
Fig2.a.wpli = plots.cf.cofiring(wpli, opt,...
'appendAttributes', {'animal', 'generateH'});
dump.Fig2 = Fig2;

% ----------------
% Fields created:
% ----------------
%.withhpc_pairs            ; % literally list of correlations between source and target=hpc
%.withpfc_pairs            ; % literally list of correlations between source and target=pfc, but separable for each batch
%.all_pairs_withhpc        ; % same as above, but all in one vector
%.all_pairs_withpfc        ; % same as above, but all in one vector
%.mean_corrwithhpc         ; % mean of all_pairs_withhpc
%.mean_corrwithpfc         ; % mean of all_pairs_withpfc
%.std_corrwithhpc          ; % std of all_pairs_withhpc
%.std_corrwithpfc          ; % std of all_pairs_withpfc
%.mean_withhpccorr_pattern ; % mean of withhpc_pairs
%.mean_withpfccorr_pattern ; % mean of withpfc_pairs
%.prophpc                  ; % properties of hpc

%% 
% The co-firing of hpc and pfc neurons during different activity patterns, with 
% average plotted

% TODO:
% 1. Add sig to figs
% 2. Redo with different normalization
% 3. Check standard coherence w/ fresh chronux
% 4. Show ranges with quantile plot

Fig2.a.all  = plots.cf.plotCofiring(Fig2.a.all, Option, 'figAppend', 'all');
Fig2.a.spec = plots.cf.plotCofiring(Fig2.a.spec,Option, 'figAppend', 'spec');
Fig2.a.coh  = plots.cf.plotCofiring(Fig2.a.coh, Option, 'figAppend', 'coh');
Fig2.a.wp   = plots.cf.plotCofiring(Fig2.a.wpli,Option, 'figAppend', 'wp');

plots.cf.plotCofiring(Fig2.a.all, Option, 'Normalization', 'probability', 'figAppend', 'all');
plots.cf.plotCofiring(Fig2.a.spec,Option, 'Normalization', 'probability', 'figAppend', 'spec');
plots.cf.plotCofiring(Fig2.a.coh, Option, 'Normalization', 'probability', 'figAppend', 'coh');
plots.cf.plotCofiring(Fig2.a.wpli,Option, 'Normalization', 'probability', 'figAppend', 'wp');

plots.cf.plotCofiring(Fig2.a.all, Option, 'Normalization', 'cdf', 'figAppend', 'all');
plots.cf.plotCofiring(Fig2.a.spec,Option, 'Normalization', 'cdf', 'figAppend', 'spec');
plots.cf.plotCofiring(Fig2.a.coh, Option, 'Normalization', 'cdf', 'figAppend', 'coh');
plots.cf.plotCofiring(Fig2.a.wpli,Option, 'Normalization', 'cdf', 'figAppend', 'wp');

%%
close all

%% 
% Example 
% the difference in cofiring between hpc-hpc pairs and hpc-pfc pairs
% are significant for all the activity patterns
tic
plots.cf.plotCfExampByDirection(Fig2.a.all,  Patterns, Option, "figAppend", 'all');
plots.cf.plotCfExampByDirection(Fig2.a.spec, Patterns, Option, "figAppend", 'spec');
plots.cf.plotCfExampByDirection(Fig2.a.coh,  Patterns, Option, "figAppend", 'coh');
plots.cf.plotCfExampByDirection(Fig2.a.wpli, Patterns, Option, "figAppend", 'wp');

plots.cf.plotCfExampByDirection(Fig2.a.all,  Patterns, Option, "figAppend", 'all',  'Normalization', 'probability');
plots.cf.plotCfExampByDirection(Fig2.a.spec, Patterns, Option, "figAppend", 'spec', 'Normalization', 'probability');
plots.cf.plotCfExampByDirection(Fig2.a.coh,  Patterns, Option, "figAppend", 'coh',  'Normalization', 'probability');
plots.cf.plotCfExampByDirection(Fig2.a.wpli, Patterns, Option, "figAppend", 'wp',   'Normalization', 'probability');
disp("Time to plot: " + toc)
%%

dump.Fig2 = Fig2;

%% 
% "These weak correlations indicate that only a small fraction of a neuron’s 
% response variability can be explained by another individual neuron"
% 
% B: Explained Variance
tic
Fig2.b = plots.pred.var.regionalexplained(Patterns, Option, ...
                        'appendAttributes', {'animal', 'generateH','name'});
Fig2.b = plots.pred.var.plotexplained(Fig2.b, Option, ...
                        'figAppend', 'all');
Fig2.b = plots.pred.var.plotexplained(Fig2.b, Option, ...
                        'figAppend', 'all-log',...
                        'yscale', 'log');
Fig2.datetime = datetime('now');
disp("Time to plot: " + toc)
%% 
close all

% Fieldss created:
% ----------------
% .r_withhpc_partitions; % firing prediction with hpc
% .r_withpfc_partitions; % firing prediction with pfc
% .single_pred_with_hpc; % single neuron firing prediction with hpc
% .single_pred_with_pfc; % single neuron firing prediction with pfc
% .patternVarExplained_hpc; % variance explained with hpc
% .patternVarExplained_pfc; % variance explained with pfc

%% Stupid test - possibly erase - not sure why this was done
% Check how dist of corr differs from dist of pred score 

% Blurt out the mean predictability of each target neuron at each condition
% (each cell holds size 1 x nNeuron)
% Stupid test - are these distributions different?
r_withhpc_patterns = [Fig2.b.r_withhpc_patterns{:}];
r_withpfc_patterns = [Fig2.b.r_withpfc_patterns{:}];
[h_hpc, p_hpc] = kstest2(r_withhpc_patterns, Fig2.a.all.all_pairs_withhpc);
[h_pfc, p_pfc] = kstest2(r_withhpc_patterns, Fig2.a.all.all_pairs_withpfc);


%% 
% Predicition Performance

% TODO: Series of heatmaps to examine corr and pred

%% GET PRED PERFORMANCE TAVLWE
[combinedPatternsTable, singlePredTable] = ...
    table.analyses.allAnimPredTable(Patterns, Option);


%% PERFORMANCE PLOTS
fig('Prediction, HPC to HPC vs HPC to PFC'); clf;
plots.pred.regionalPerf(combinedPatternsTable, singlePredTable)
% Run again for each genH
for genh = unique(combinedPatternsTable.genH)'
    fig(sprintf('Prediction, HPC to HPC vs HPC to PFC - %s', genh{1})); clf;
    plots.pred.regionalPerf(combinedPatternsTable, ...
        singlePredTable, 'select_genH', genh)
end

% Individual samples
plots.pred.plotIndivSamples(combinedPatternsTable, Option);
plots.pred.plotIndivSamples(combinedPatternsTable, Option, 'power')
plots.pred.plotIndivSamples(combinedPatternsTable, Option, 'coherence')
plots.pred.plotIndivSamples(combinedPatternsTable, Option, 'wpli')
%%

% print stats
formatSpec1 = "%s: %0.3f±%0.3f";
disp("Prediction Performance on Average")
sprintf(formatSpec1,Patterns(1,1).directionality,mean_withhpc,std_withhpc)
sprintf(formatSpec1,Patterns(2,1).directionality,mean_withpfc,std_withpfc)
if useSinglePrediction
    disp('single source prediction median')
    formatSpec2 = "%s: %0.5e";
    sprintf(formatSpec2,Patterns(1,1).directionality,median_singlehh)
    sprintf(formatSpec2,Patterns(1,2).directionality,median_singlehp)
end

%% Figures --- Prediction Performance with PredDims
plots.pred.withPredDims
% See python script for the rest of the figures Notebooks/python/dimPred.py

%% (Factor analysis)
% histogram of all predicition perf
hpc_theta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-hpc" & T.patternType == "theta";
hpc_delta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-hpc" & T.patternType == "delta";
hpc_ripple = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-hpc" & T.patternType == "ripple";


pfc_theta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.patternType == "theta";
pfc_delta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.patternType == "delta";
pfc_ripple = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.patternType == "ripple";

nPartitions = max([Patterns.iPartition]); 
nPatternsAndControls = Option.nPatternAndControl;
temp1 = zeros(3,nPartitions);
temp2 = zeros(3,nPartitions);
temp3 = zeros(3,nPartitions);

if isfield(Patterns, "factorAnalysis")
    for j = 1:nPartitions
        for k = 1:nPatternsAndControls
            temp2(k,j) = Patterns(j,hpc,k).factorAnalysis.qOpt;
            temp3(k,j) = Patterns(j,pfc,k).factorAnalysis.qOpt;
            temp1(k,:) = temp2(k,:)./temp3(k,:);
        end
    end
end
ratios      = mean(temp1,2);
hpcQoptDims = mean(temp2,2);
pfcQoptDims = mean(temp3,2);


%% Figure 6 (Figure 4 in Thesis)
%% 
% *Pattern Specific* 
plots.pred.removePred
writetable(rt, figuredefine("tables", "rt_zscr="+zscr+".csv"), 'WriteRowNames',true)

%  Remove tuples of (targetArea,pattern) from a given (targetArea, pattern)
disp("Running dimension removal with zscore=" + zscr);
% plot with grm
if zscr; figAppend = "zscore"; else; figAppend = ""; end
plots.grm.dimensionRemoved(rt, "spec", 'figAppend', figAppend);
plots.grm.dimensionRemoved(rt, "coh", 'figAppend', figAppend);
plots.grm.dimensionRemoved(rt, "wpli", 'figAppend', figAppend);
plots.subspace.pred_dim_rem.oneDimesionRemovedClustermap; % changed this:: ask ziyi if she has missing function
% plots.subspace.pred_dim_rem.plotDimensionRemoval_perPattern; # not used
plots.subspace.pred_dim_rem.plotDimensionRemoval_perPatternbyDirection; % important

%% 
% Now let's see how similar those curves are in dimension reduced space


%% Figure 7 (REQUIRES FACTOR ANALYSIS)
% 
% A: Dominant dimension predicitions

% numUsedForPrediction = min(nTarget,nSource);
% % make the averaged version
% curr_cvLoss = cell(10,2,3);
% curr_qOptDim = cell(10,2,3);
% for p = 1:10
%     for i = 1:nPatterns
%         for j = 1:2
%             if ~Patterns(p,j,i).singularWarning
%                 curr_cvLoss{p,j,i} = Patterns(p,j,i).factorAnalysis.cvLoss;
%                 curr_qOptDim{p,j,i} = Patterns(p,j,i).factorAnalysis.optDimFactorRegress;
%             end
%         end
%     end
% end
%%
try
    fig("Model performance versus factor analysis dimension")
    clf
    full_model_performance = [];

    for i = 1:nPatterns
        for j = 1:2
            subplot(3,2,2*(i-1)+j)
            
            full_model = plots.plotPredictiveDimensions(numUsedForPrediction,curr_cvLoss(:,j,i),'optDim',curr_qOptDim, ...
                "mode", "fa");
            
            full_model_performance = [full_model_performance, full_model];
            xlim([0,12.5])
            hold on
            
            plot(1, full_model,'^');
            ylim([0, full_model+0.1])
            if j == 1
                ax1 = gca;
            else
                ax2 = gca;
            end
            title([Patterns(p,j,i).name Patterns(p,j,i).directionality])
        end
        linkaxes([ax1,ax2],'y')
        
    end
catch
    warning("Factor analysis not run")
end
