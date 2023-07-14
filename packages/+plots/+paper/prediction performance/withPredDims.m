%% this script will calculate and plot the preformance as predicitive dimensions
% increases

%%
% make the averaged version
curr_cvLoss = cell(nAnimal, nMethods, nAnimalPartition, 2, nPatterns);
curr_rrDim = cell(nAnimal, nMethods, nAnimalPartition, 2, nPatterns);
for a = 1:nAnimal
for m = 1:nMethods
    for p = 1:nAnimal
        for i = 1:nPatterns
            for j = 1:2
                curr_cvLoss{a,m,p,j,i} = Patterns(a,m,p,j,i).rankRegress.cvLoss;
                curr_rrDim{a,m,p,j,i} = Patterns(a,m,p,j,i).rankRegress.optDimReducedRankRegress;
            end
        end
    end
    
end
end

%% calculate/plot
figure(701)
clf
mean_full_model_performance = [];

for m = 1:1
    for i = 1:nPatterns
        for j = 1:2
            subplot(1,3,i)
            
            if j == 1
                full_model = plotPredictiveDimensions(numDimsUsedForPrediction{2},...
                    curr_cvLoss(2,m,1,j,i), "optDim", curr_rrDim(2,m,1,j,i),"mode", "rr", "color","blue", "normalized", true);
            else
                full_model = plotPredictiveDimensions(numDimsUsedForPrediction{2},...
                    curr_cvLoss(2,m,1,j,i), "optDim", curr_rrDim(2,m,1,j,i),"mode", "rr" ,"normalized", true);
            end
            hold on
            mean_full_model_performance = [mean_full_model_performance, full_model];
            xlim([0,8])

            if j == 1
                ax1 = gca;
            else
                ax2 = gca;
            end
            title([Patterns_AllAnimals(m,1,j,i).name ])
            xlim([0,10])
        end
        linkaxes([ax1,ax2],'y')
        
    end
    hold on
end

% how to mark legends?
% how to do different animals?
legend ("hpc-hpc", "hpc-pfc")