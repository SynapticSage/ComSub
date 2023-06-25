function [Patterns, Patterns_overall] = partition(r, Option) 
% partition - Partition data into source and target areas
% 
% Parameters:
%   r - A struct containing raw and processed data
%   Option - A struct containing options for the partitioning
% 
% Returns:
%   Patterns - A struct containing the partitioned data (partitions of
%               neurons; subsampled sets of them)
%   Patterns_overall - A struct WITHOUT parititioning neurons
%
% ----------------------------
% Note:
% ----------------------------
% when the partition is three-ways, direction==1 means same target/source pair
% and direction==2 means diff target/source pair
% ----------------------------
disp('Partitioning data...')
tic
 
% Initialize the scaffold of the pattern struct
Patterns         = initPatternStruct();
Patterns_overall = initPatternStruct();

% ------------------------------
% Place paritioned data properly
% ------------------------------
for iPartition = 1:Option.numPartition

    % Split cells into source and target areas
    % [s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index,...
    %  s_pfc_index, t_hpc_index, t_pfc_index, directionality] = ...
    %     trialSpikes.split.legacyRun(r, Option);
    split_info = trialSpikes.split.run(r, Option);

    for i = 1:numel(r.windowInfo.cellOfWindows)
        for j = 1:numel(split_info.directionality)

            
            % Parse directionality
            directionality = split_info.directionality(j);
            sourcetarg = directionality.split('-');
            source = sourcetarg(1);
            target = sourcetarg(1);
            s_dat = split_info.source;
            s_ind = split_info.source_index;
            t_dat = split_info.target{j, 1, i};
            t_ind = split_info.target_index(j, :);
            
            % Assign x_source and x_target
            Patterns(iPartition,j,i).X_source  = s_dat;
            Patterns(iPartition,j,i).X_target  = t_dat;
            Patterns(iPartition,j,i).X_time    = reshape(r.trialTimes{i}',1,[]);
            % Assign index_source and index_target
            Patterns(iPartition,j,i).index_source = s_ind;
            Patterns(iPartition,j,i).index_target = t_ind;
            % Assign directionality
            Patterns(iPartition,j,i).directionality = directionality;
            % Assign pattern name
            try
                Patterns(iPartition,j,i).name = Option.patternNames(i);
            catch 
                Patterns(iPartition,j,i).name = "Pattern " + i;
            end

            if iPartition == 1
                if source == "hpc"
                    s_all = r.hpc.X;
                    s_ind_all = 1:size(r.hpc.X, 2);
                elseif source == "pfc"
                    s_all = r.pfc.X;
                    s_ind_all = 1:size(r.pfc.X, 2);
                else
                    error('Source area not recognized')
                end
                if target == "hpc"
                    t_all = r.hpc.X;
                    t_ind_all = 1:size(r.hpc.X, 2);
                elseif target == "pfc"
                    t_all = r.pfc.X;
                    t_ind_all = 1:size(r.pfc.X, 2);
                else
                    error('Target area not recognized')
                end
                % Assign x_source and x_target
                Patterns_overall(j,i).X_source = s_all;
                Patterns_overall(j,i).X_target = t_all;
                % Assign index_source and index_target
                Patterns_overall(j,i).index_source = s_ind_all;
                Patterns_overall(j,i).index_target = t_ind_all;
                % Assign directionality
                Patterns_overall(j,i).directionality = directionality;
                % Assign pattern name
                try
                    Patterns_overall(iPartition,j,i).name = Option.patternNames(i);
                catch 
                    Patterns_overall(iPartition,j,i).name = "Pattern " + i;
                end
            end

        end
    end
end

szCellOfWindows = squeeze(size(r.windowInfo.cellOfWindows));
if numel(szCellOfWindows) == 2 && szCellOfWindows(1) == 1
    szCellOfWindows = szCellOfWindows(2);
end
Patterns = reshape(Patterns, [Option.numPartition, numel(split_info.directionality), size(r.windowInfo.cellOfWindows)]);

disp(['Partitioning data took ', num2str(toc), ' seconds.'])
