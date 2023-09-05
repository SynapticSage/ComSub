
%%
tmp = nd.fieldGet(Patterns,'rankRegress');
if iscell(tmp)
    tmp = ndb.toNd(tmp);
end
optDim = nd.fieldGet(tmp,'optDimReducedRankRegress');
optDim = floor(median(optDim, 1));
optDim = max(optDim,[],'all');
patternnames = Option(1).patternNames(1:size(Patterns,ndims(Patterns)));

ugenH = shortcut.genH_shorter(unique([Patterns.generateH]));

%%
% Run all of them
rt = table();
for genh = ugenH
    [Pat, opt] = munge.getH(Patterns, Option, genh);
    assert(~isempty(Pat), 'No patterns found for %s', genh);
    Patterns_AllAnimals = squeeze(Pat);
    Patterns_AllAnimals = reshape(Patterns_AllAnimals, [size(Patterns_AllAnimals,1)*size(Patterns_AllAnimals,2), size(Patterns_AllAnimals,3), size(Patterns_AllAnimals,4)]);
    for iPartition = progress(1:size(Patterns_AllAnimals,2),'Title','iPartition')
        for baseDirection = 1:2
            for removeDirection = 1:2
                for basePattern = 1:numel(patternnames)
                    for removePattern = 1:numel(patternnames)
                        
                        % Base patterns
                        X_source = Patterns_AllAnimals(iPartition, baseDirection, basePattern).X_source;
                        X_target = Patterns_AllAnimals(iPartition, baseDirection, basePattern).X_target;
                        cvLoss = Patterns_AllAnimals(iPartition,   baseDirection, basePattern).rankRegress.cvLoss;
                        optDim = Patterns_AllAnimals(iPartition,   baseDirection, basePattern).rankRegress.optDimReducedRankRegress;
                        numDimsUsedForPrediction = 1:size(cvLoss,2);
                        nTarget = size(X_target,1);
                        % Remove patterns
                        B_ = Patterns_AllAnimals(iPartition, removeDirection, removePattern).rankRegress.B_;
                        % Run dim removal
                        [performance,fullmodel] = plots.sequentialRemovePredDims(X_source, X_target, B_, optDim,...
                            cvLoss, numDimsUsedForPrediction, false, "normalized", false);

                        % fractionaldim = (1:size(cvLoss,2))./double(optDim);
                        % fractionaldim = fractionaldim(:);
                        
                        % size(genh), size(iPartition), size(baseDirection), size(removeDirection), size(patternnames(basePattern)), size(patternnames(removePattern)), size(fullmodel)
                        
                        t = table(genh, iPartition, baseDirection,...
                        removeDirection, patternnames(basePattern), ...
                        patternnames(removePattern),fullmodel, ...fractionaldim,...
                        'VariableNames',{'method','partition',...
                        'baseDirection','removeDirection','basePattern',...
                        'removePattern','fullmodel', ...
                        ...'fractonaldim'...
                        });
                        t = repmat(t,numel(performance),1);
                        performance = performance(:);
                        t.dimensionRemoved = (0:numel(performance)-1)';
                        t.performance = performance(:);
                        rt = [rt; t];
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
