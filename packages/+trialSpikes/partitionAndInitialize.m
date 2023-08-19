function [Patterns, Patterns_overall] = partition(Spk, Option, varargin)
% partition - Partition data into source and target areas
% 
% Parameters:
%   Spk - A struct containing raw and processed data
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

ip = inputParser;
ip.addParameter('test', false, @islogical); % test mode, use fewer partitions
ip.parse(varargin{:});
Opt = ip.Results;

disp('Partitioning data...')
tic
const = option.constants();
 
% Initialize the scaffold of the pattern struct
Patterns         = initPatternStruct();
Patterns_overall = initPatternStruct();

% ------------------------------
% Place paritioned data properly
% ------------------------------
if Opt.test % test mode, use fewer partitions
    Option.numPartition = 2;
end
for iPartition = progress(1:Option.numPartition, 'Title', 'Partitioning data')
% for iPartition = progress(1:2, 'Title', 'Partitioning data')

    % Split cells into source and target areas
    % [s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index,...
    %  s_pfc_index, t_hpc_index, t_pfc_index, directionality] = ...
    %     trialSpikes.split.legacyRun(Spk, Option);
    split_info = trialSpikes.split.run(Spk, Option);

    for i = 1:Option.nPatternAndControl
        for j = 1:numel(split_info.directionality)

            % Parse directionality
            directionality = split_info.directionality(j);
            sourcetarg = directionality.split('-');
            source = sourcetarg(1);
            target = sourcetarg(2);
            s_dat = split_info.source{i};
            s_ind = split_info.source_index;
            t_dat = split_info.target{j, 1, i};
            t_ind = split_info.target_index(j, :);
            % Do a "checksum"-ish move if source==target
            if source == target
                assert(numel(intersect(s_ind, t_ind)) == 0, ...
                    'If source==target, then source and target indices should be disjoint')
            end
            
            % Assign x_source and x_target
            assert(~clean.isRankDeficient(s_dat), ...
                'Source data is rank deficient');
            Patterns(iPartition,j,i).X_source  = s_dat;
            Patterns(iPartition,j,i).X_target  = t_dat;
            assert(all(size(s_dat) > 0), 'Source data is empty');
            assert(all(size(t_dat) > 0), 'Source data is empty');
            Patterns(iPartition,j,i).X_time    = reshape(Spk.trialTimes{i}',1,[]);
            % Assign index_source and index_target
            Patterns(iPartition,j,i).index_source = s_ind;
            Patterns(iPartition,j,i).index_target = t_ind;
            ind_ogsource = Spk.hpc.to_orig(s_ind);
            Patterns(iPartition,j,i).index_og_source = ind_ogsource;
            if directionality == "hpc-hpc"
                ind_ogtarget = Spk.hpc.to_orig(t_ind);
                Patterns(iPartition,j,i).index_og_target = ind_ogtarget;
            elseif directionality == "hpc-pfc"
                ind_ogtarget = Spk.pfc.to_orig(t_ind);
                Patterns(iPartition,j,i).index_og_target = ind_ogtarget;
            else
                error('Not implemented')
            end
            % Patterns(iPartition,j,i).muFR_X_source = Spk.muFR(ind_ogsource);
            % Patterns(iPartition,j,i).muFR_X_target = Spk.muFR(ind_ogtarget);
            % Patterns(iPartition,j,i).stdFR_X_source = Spk.stdFR(ind_ogsource);
            % Patterns(iPartition,j,i).stdFR_X_target = Spk.stdFR(ind_ogtarget);
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
                    s_all = Spk.hpc.X{i};
                    s_ind_all = 1:size(s_all, 1);
                elseif source == "pfc"
                    s_all = Spk.pfc.X{i};
                    s_ind_all = 1:size(s_all, 1);
                else
                    error('Source area not recognized')
                end
                if target == "hpc"
                    t_all = Spk.hpc.X{i};
                    t_ind_all = 1:size(t_all, 1);
                elseif target == "pfc"
                    t_all = Spk.pfc.X{i};
                    t_ind_all = 1:size(t_all, 1);
                else
                    error('Target area not recognized')
                end
                if source ~= target
                    assert(numel(s_ind_all) ~= numel(t_ind_all), 'Source and target areas are the same')
                end
                directionality = source + "-" + target;
                % Assign x_source and x_target
                assert(~clean.isRankDeficient(s_all));
                Patterns_overall(j,i).X_source = s_all;
                Patterns_overall(j,i).X_target = t_all;
                Patterns_overall(j,i).X_time = reshape(Spk.trialTimes{i}',1,[]);
                % Assign index_source and index_target
                Patterns_overall(j,i).index_source = s_ind_all;
                Patterns_overall(j,i).index_target = t_ind_all;
                % Assign directionality
                Patterns_overall(j,i).directionality = directionality;
                % Assign pattern name
                try
                    Patterns_overall(j,i).name = Option.patternNames(i);
                    % Full name = PatternNamesFull + directionality
                    Patterns_overall(j,i).nameFull = Option.patternNamesFull(i) + ...
                        "-(" + directionality + ")";
                catch 
                    Patterns_overall(j,i).name = "Pattern " + i;
                    Patterns_overall(j,i).nameFull = "Pattern " + i + ...
                        "-(" + directionality + ")";
                end
            end

        end
    end
end

% Add overall for ALL TIME
last = Option.nPatternAndControl + 1;
areaPerNeuron     = Spk.areaPerNeuron;
sessionTypePerBin = Spk.sessionTypePerBin;
for j = 1:numel(split_info.directionality)
    if j  == const.HPC
        s_all = Spk.spikeRateMatrix(areaPerNeuron == "CA1", ... 
            sessionTypePerBin == 1);
        t_all = s_all;
        s_inds = Spk.celllookup(areaPerNeuron == "CA1",:).index_by_region;
        t_inds = s_inds;
        source = "hpc";
        target = "hpc";
    else
        s_all = Spk.spikeRateMatrix(areaPerNeuron == "CA1", ... 
            sessionTypePerBin == 1);
        t_all = Spk.spikeRateMatrix(areaPerNeuron == "PFC", ... 
            sessionTypePerBin == 1);
        s_inds = Spk.celllookup(areaPerNeuron == "CA1",:).index_by_region;
        t_inds = Spk.celllookup(areaPerNeuron == "PFC",:).index_by_region;
        source = "hpc";
        target = "pfc";
    end
    % Assign x_source and x_target
    Patterns_overall(j, last).X_source = s_all;
    Patterns_overall(j, last).X_target = t_all;
    % Assign index_source and index_target
    Patterns_overall(j, last).index_source = s_inds;
    Patterns_overall(j, last).index_target = t_inds;
    % Assign directionality
    Patterns_overall(j, last).directionality = directionality;
    % Assign pattern name
    Patterns_overall(j, last).name = "Overall";
    time = Spk.timeBinMidPoints(sessionTypePerBin == 1);
    Patterns_overall(j, last).X_time = reshape(time',1,[]);
    Patterns_overall(j, last).source = source;
    Patterns_overall(j, last).target = target;
end

szCellOfWindows = squeeze(size(Spk.spikeSampleMatrix));
if numel(szCellOfWindows) == 2  && szCellOfWindows(1) == 1
    szCellOfWindows = szCellOfWindows(2);
end

Patterns = reshape(Patterns, ...
    [Option.numPartition, numel(split_info.directionality), ...
        size(Spk.spikeSampleMatrix)]);

disp(['Partitioning data took ', num2str(toc), ' seconds.'])
