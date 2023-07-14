% this script calculates the single cell firing predicition of all patterns
% across methods
%% initialize data sturctures
single_pred_with_hpc = cell(nMethods, nPatterns, nAnimalPartition);
single_pred_with_pfc = cell(nMethods, nPatterns, nAnimalPartition);

%% loop
for m = 1:nMethods
    for i = 1:nPatterns
        for p = 1:nAnimalPartition
            
            
            curr_source = Patterns(m,p,1,i).X_source';
            nSource = size(curr_source,2);
            for d = 1:2
                curr_target = Patterns(m,p,d,i).X_target';
                for k = 1:nSource
                    curr_singleB = Patterns(m,p,d,i).rankRegress.singlesource_B{k};
                    if ~isempty(curr_singleB)
                        curr_singlesource = curr_source(:,k);
                        [singlepattern,~] = calculateVarianceExplained(curr_singlesource, curr_target, curr_singleB);
                        if d == 1
                            [single_pred_with_hpc] = [single_pred_with_hpc singlepattern];
                        else
                            [single_pred_with_pfc] = [single_pred_with_pfc singlepattern];
                        end
                    end
                end
            end
        end
    end
end
median_singlehh = median(single_pred_with_hpc(~isnan(single_pred_with_hpc)));
temp1 = single_pred_with_pfc(~isinf(single_pred_with_pfc));
temp1 = temp1(~isnan(temp1));
median_singlehp = median(temp1);