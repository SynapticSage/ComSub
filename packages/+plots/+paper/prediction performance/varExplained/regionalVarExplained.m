% this script calculates the regional firing predicition of all patterns
% across methods
%% initialize data sturctures
r_withhpc_patterns = cell(nMethods, nPatterns, nAnimalPartition);
r_withpfc_patterns = cell(nMethods, nPatterns, nAnimalPartition);

single_pred_with_hpc = cell(nMethods, nPatterns, nAnimalPartition);
single_pred_with_pfc = cell(nMethods, nPatterns, nAnimalPartition);

%% loop
for m = 1:nMethods
    for i = 1:nPatterns
        for p = 1:nAnimalPartition
            animalNo = floor(p/(nPartition+1)+1);
            nCurrSource = nSource(animalNo);
            nCurrTarget = nTarget(animalNo);
            single_pred_with_hpc{m}{i}{p} = [];
            single_pred_with_pfc{m}{i}{p} = [];

            r_withhpc_patterns{m}{i}{p} = zeros(1,nCurrTarget);
            r_withpfc_patterns{m}{i}{p} = zeros(1,nCurrTarget);
            
            curr_source = (Patterns_AllAnimals(m,p,1,i).X_source)';
            curr_targethpc = (Patterns_AllAnimals(m,p,hpc,i).X_target)';
            curr_targetpfc = (Patterns_AllAnimals(m,p,pfc,i).X_target)';
            
            curr_B_hpc = Patterns_AllAnimals(m,p,hpc,i).rankRegress.B;
            curr_B_pfc = Patterns_AllAnimals(m,p,pfc,i).rankRegress.B;
            
            [patternhpc, meanhpc, ~] = calculateVarianceExplained...
                                      (curr_source, curr_targethpc, curr_B_hpc);
            [patternpfc, meanpfc, ~] = calculateVarianceExplained...
                                      (curr_source, curr_targetpfc, curr_B_pfc);
            
            patternVarExplained_hpc(m,i,p) = meanhpc;
            patternVarExplained_pfc(m,i,p) = meanpfc;
            
            r_withpfc_patterns{m}{i}{p} = patternpfc;
         
            r_withhpc_patterns{m}{i}{p} = patternhpc;
            
            
%             for k = 1:nCurrSource
%                 curr_singleB = Patterns_AllAnimals...
%                                (m,p,1,i).rankRegress.singlesource_B{k};
%                 if ~isempty(curr_singleB)
%                      curr_singlesource = curr_source(:,k);
%                     [single_hpc,~] = calculatePredictionPerformance...
%                                               (curr_singlesource, curr_targethpc, curr_singleB);
%                     [single_pfc,~] = calculatePredictionPerformance...
%                                               (curr_singlesource, curr_targetpfc, curr_singleB);
%                     single_pred_with_hpc{m}{i}{p} = ...
%                                      [single_pred_with_hpc{m}{i}{p}, single_hpc];
%                     single_pred_with_pfc{m}{i}{p} =...
%                                      [single_pred_with_pfc{m}{i}{p}, single_pfc];
%                 end
%             end
        end
    end
end
