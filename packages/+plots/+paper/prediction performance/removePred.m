
%%
tmp = nd.fieldGet(Patterns,'rankRegress');
optDim = nd.fieldGet(tmp,'optDimReducedRankRegress');
optDim = floor(median(optDim, 1));
optDim = max(optDim,[],'all');

%%
% Run all of them
rt = table();
for iMethod = 1:1
    for iPartition = progress(1:size(Patterns_AllAnimals,2),'Title','iPartition')
        for baseDirection = 1:size(Patterns_AllAnimals,3)
            for removeDirection = 1:size(Patterns_AllAnimals,3)
                for basePattern = 1:numel(patternnames)
                    for removePattern = 1:numel(patternnames)
                        
                        % Base patterns
                        X_source = Patterns_AllAnimals(iMethod, iPartition, baseDirection, basePattern).X_source;
                        X_target = Patterns_AllAnimals(iMethod, iPartition, baseDirection, basePattern).X_target;
                        cvLoss = Patterns_AllAnimals(iMethod, iPartition,   baseDirection, basePattern).rankRegress.cvLoss;
                        numDimsUsedForPrediction = 1:size(cvLoss,2);
                        % Remove patterns
                        B_ = Patterns_AllAnimals(iMethod, iPartition, removeDirection, removePattern).rankRegress.B_;
                        % Run dim removal
                        [performance,fullmodel] = plots.sequentialRemovePredDims(X_source, X_target, B_, optDim,...
                            cvLoss, numDimsUsedForPrediction, false, "normalized", false);
                        
                        
                        T = table(iMethod, iPartition, baseDirection, removeDirection, patternnames(basePattern), patternnames(removePattern),fullmodel,...
                            'VariableNames',{'method','partition','baseDirection','removeDirection','basePattern','removePattern','fullmodel'});
                        T = repmat(T,numel(performance),1);
                        T.dimensionRemoved = (0:numel(performance(:))-1)';
                        T.performance = performance(:);
                        rt = [rt; T];
                    end
                end
            end
        end
    end
end


%%
rt.sameDirection = rt.baseDirection == rt.removeDirection;
dirlabel = ["remove other area","remove same area"]';
rt.sameDirectionLabel = dirlabel(rt.sameDirection+1);
dirlabel = ["hpc-hpc","hpc-pfc"]';
rt.targetArea = dirlabel(rt.baseDirection);
rt.removeTargetArea = dirlabel(rt.removeDirection);
rt.basePatternLabel = rt.basePattern.replace('theta','$\theta$').replace('delta','$\delta$').replace('ripple','SPW-R');
rt.removePatternLabel = rt.removePattern.replace('theta','$\theta$').replace('delta','$\delta$').replace('ripple','SPW-R');