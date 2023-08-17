% Tries to display raw data like the 2016 figurevg
% Inputs:
%   Spk - struct of spike data
%   Events - struct of event data
%   Patterns_overall - struct of CCA data
%   behavior - table of behavior data
%
% Notes
% - Quantile filter by different behavior var 
% - occ norm  
% - cca, dim 3 
% - color cells by component 
% - actual efizz 🤷 
% - window eeg 
%    - color by component 
% - subplot grid with behavior for selected period (1 traj per subplot)  
close all

%% SETTINGS
delta_band = [0.5, 4];
theta_band = [6, 12];
ripple_band = [150, 250];
gap_thresh = 0.5; % Adjust as necessary
use_rescale_avg = true;
sortprop = 'dirlindist';
sortdir = 'descend';
dim3 = [];
colorbycomp = true;
wins = [1, 3]; % which windows to plot
shadeOption = true; % show windows : set to true to draw rectangles, and false to skip
lfp_fields = ["theta", "ripple", "data"];
use_fft_avg = true; % otherwise use lfp
% For behavior based plotting
sets_wanna_plot = {...
["u",1,"X_time"], ["v",1,"X_time"], ["u",2,"X_time"], ["v",2,"X_time"], ["u",3,"X_time"], ["v",3,"X_time"],... 
...["theta_hpc",1,"t"], ["theta_pfc",1,"t"], ["ripple_hpc",1,"t"],["ripple_pfc",1,"t"],...
...["theta_wpli_hpc",1,"t"], ["theta_wpli_pfc",1,"t"], ["ripple_wpli_hpc",1,"t"],["ripple_wpli_pfc",1,"t"],...
...["delta_wpli_hpc",1,"t"], ["delta_wpli_pfc",1,"t"], ["theta_wpli_hpc",2,"t"], ["theta_wpli_pfc",2,"t"],...
};
const = option.constants();
all_animals = const.all_animals;
all_animals = ["ZT2" setdiff(all_animals, "ZT2")];
all_animals = setdiff(all_animals, ["ZT2", "ER1","JS13"]);
compcolors = [
    0,   128, 128; ...  % Strong Teal
    255, 165, 0;   ...  % Strong Orange
    139, 0,   139; ...  % Contrast Color (Purple/Violet)
    128, 200, 200; ...  % Light Teal
    255, 218, 176; ...  % Light Orange
    214, 161, 214; ...  % Light Purple/Violet
    64,  152, 152; ...  % Mid Teal
    255, 194, 120; ...  % Mid Orange
    177, 100, 177  ...  % Mid Purple/Violet
] / 255;  % Dividing by 255 to convert values into MATLAB's RGB range [0,1]
plots.plot_discrete_colors(compcolors);
sgtitle("Window components");


%% LOAD DATA
for animal = all_animals
matfile_name=@(animal) hashdefine(animal + "_fromCoherencefromRipTimes_zscore_true_midpattern_true_mostRecentState.mat");
checkpoint = matfile_name(animal);

if ~exist("efizz", "var") && ~exist('Option') || animal ~= Option.animal
    commsubspaceToPath
    clear efizz Spk Option Patterns_overall behavior
    Option = store.load_checkpoint(checkpoint, 'Option');
    Option.save=false;
    Option.analysis.rankRegress=false;
    Option.analysis.cca=false;
    TheScript;
    figuredefine("-clearpermfolder")
    load(Option.animal + "spectralBehavior.mat");
    load(Option.animal + "avgeeg.mat");
    avgeeg = ndb.toNd(avgeeg);
    lfp.hpc = munge.oneEEGstruct(avgeeg(:,:,1)); % TODO: how does this handle areas?
    lfp.pfc = munge.oneEEGstruct(avgeeg(:,:,2));
    lfp.time = mean([lfp.hpc.time; lfp.pfc.time], 1);
    clear avgeeg
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    behavior.dirlindist = sign(behavior.trajbound-0.5) .* behavior.lindist;
    Patterns_overall = analysis.cca(Patterns_overall, Option);
end
pattern_overall_ind = numel(Patterns_overall);
savefolder_opt = "_showwin=" + num2str(shadeOption) + "_colorbycomp=" + num2str(colorbycomp) + "_dim3=" + num2str(dim3) + "_sortprop=" + sortprop + "_sortdir=" + sortdir + "_use_fft_avg=" + use_fft_avg;
if Option.midpattern
    savefolder = figuredefine("mipattern=true","2016figure", savefolder_opt);
else
    savefolder = figuredefine("2016figure",savefolder_opt);
end
if ~exist(fileparts(savefolder), 'dir')
    mkdir(savefolder);
end
if ~exist(savefolder, 'dir')
    mkdir(savefolder);
end
%% Specify the epoch, traj, trajall, trajclass, trajbound (use [] if not specifying)
for epoch_option     = progress(2:2:16,'Title','Epoch'); % e.g., 'EpochName' %FOR_EPOCH
% for trajbound_option = [0 1]; % 0 or 1
for trajbound_option = 0:1; % 0 or 1 FOR_TRAJBOUND
traj_option      = []; % a specific trajectory index within an epoch
trajall_option   = []; % a specific trajectory index over all epochs
trajclass_option = []; % 1 through 4
%% -- INDEXING PROPERTIES --

% Find frequency indices
delta_idx  = find(efizz.f >= delta_band(1) & efizz.f <= delta_band(2));
theta_idx  = find(efizz.f >= theta_band(1) & efizz.f <= theta_band(2));
ripple_idx = find(efizz.f >= ripple_band(1) & efizz.f <= ripple_band(2));

% Compute average power in the theta and ripple bands for both regions
avg.theta_hpc   = mean(efizz.S1(:, theta_idx), 2);
avg.ripple_hpc  = mean(efizz.S1(:, ripple_idx), 2);
avg.theta_pfc   = mean(efizz.S2(:, theta_idx), 2);
avg.ripple_pfc  = mean(efizz.S2(:, ripple_idx), 2);
avg.theta_cavg  = mean(efizz.Cavg(:, theta_idx), 2);
avg.ripple_cavg = mean(efizz.Cavg(:, ripple_idx), 2);
avg.delta_cavg  = mean(efizz.Cavg(:, delta_idx), 2);
avg.theta_wpli  = mean(efizz.wpli_avg(:, theta_idx), 2);
avg.ripple_wpli = mean(efizz.wpli_avg(:, ripple_idx), 2);
avg.delta_wpli  = mean(efizz.wpli_avg(:, delta_idx), 2);

% Define hpc cells
cells.hpc = Spk.areaPerNeuron  == "CA1";
cells.pfc = Spk.areaPerNeuron  == "PFC";

%% TIME :: Finding Valid Gaps in Data


% Filter the behavior table based on the provided options
selected_rows = behavior;
if ~isempty(epoch_option)
    selected_rows = selected_rows(selected_rows.epoch == epoch_option, :);
end
if ~isempty(traj_option)
    selected_rows = selected_rows(selected_rows.traj == traj_option, :);
end
if ~isempty(trajall_option)
    selected_rows = selected_rows(selected_rows.trajall == trajall_option, :);
end
if ~isempty(trajclass_option)
    selected_rows = selected_rows(selected_rows.trajclass == trajclass_option, :);
end
if ~isempty(trajbound_option)
    selected_rows = selected_rows(selected_rows.trajbound == trajbound_option, :);
end
if isempty(selected_rows)
    warning("No data for " + animal + " animal " + epoch_option + " epoch " + trajbound_option + " trajbound");
    continue
end

% Extract start and stop times from the filtered rows
start_time = min(selected_rows.time);
stop_time  = max(selected_rows.time);
% Find time gaps in behavior's time field
gaps = find(diff(selected_rows.time) > gap_thresh);
% Use gaps to derive ranges (begin and end of each segment without a gap)
g = [selected_rows.time(gaps)'; selected_rows.time(gaps+1)'];
initial = [start_time; g(1,1)];
final = [g(2,end); stop_time];
G = [g(2,1:end-1); g(1,2:end)];
ranges = [initial, G, final];
total_time = sum(ranges(2,:) - ranges(1,:));
disp("Total time: " + total_time + " seconds");
% Calculate the actual time gap between each range and the previous range's end
% time_gaps = diff(ranges(2,:)) - diff(ranges(1,:));
% time_gaps = ranges(2,1:end-1) - ranges(1,2:end);
time_gaps = ranges(1,2:end) - ranges(2,1:end-1);
% Calculate the cumulative sum of these gaps, which will be the offset for each range
cumulative_gaps = [0, cumsum(time_gaps)];
Events.wincenter = cellfun(@(x)mean(x, 2), Events.cellOfWindows, 'UniformOutput', false);

% figure(fig("Time segments"));clf;
% tiledlayout('flow'); nexttile
% plot(1:size(ranges, 2), ranges(1,:), 'bo', 'DisplayName', 'Start Times');
% yline(ranges(1,:), ':', 'Color', 'b');
% hold on;
% plot(1:size(ranges, 2), ranges(2,:), 'rx', 'DisplayName', 'End Times');
% yline(ranges(2,:), ':', 'Color', 'r');
% xlabel('Time');
% ylabel('Segment Number');
% % legend();
% indsof = @(x) 1:numel(x);
% hold on; plot(Spk.timeBinMidPoints(ind.spike), 'kx', 'DisplayName', 'Spikes');
% hold on; plot(selected.spike.timeBinMidPoints, 'gx', 'DisplayName', 'Selected Spikes');
% nexttile
% plot(Spk.timeBinMidPoints(ind.spike), selected.spike.timeBinMidPoints, 'kx');

b=behavior(behavior.time > ranges(1,1) & behavior.time < ranges(2,1), :);
b.trajbound
b.time

%% TIME :: Indexing
% Initialize struct for indices
ind = struct();
% Extract indices for each segment and data type
for i = 1:size(ranges, 2)
    % For Spk
    ind.spike{i} = find(Spk.timeBinMidPoints >= ranges(1, i) & Spk.timeBinMidPoints <= ranges(2, i));
    % For efizz
    ind.efizz{i} = find(efizz.t >= ranges(1, i) & efizz.t <= ranges(2, i));
    % For Patterns_overall
    ind.pattern{i} = find(Patterns_overall(end).X_time >= ranges(1, i) & Patterns_overall(end).X_time <= ranges(2, i));
    ind.lfp{i} = find(lfp.time >= ranges(1, i) & lfp.time <= ranges(2, i));
    for j = 1:numel(Events.wincenter)
        ind.wins{j}{i} = find(Events.wincenter{j} >= ranges(1, i) & Events.wincenter{j} <= ranges(2, i));
    end
end
ind.spike   = cat(2, ind.spike{:});
ind.efizz   = cat(2, ind.efizz{:});
ind.pattern = cat(2, ind.pattern{:});
ind.lfp = cat(2, ind.lfp{:});
for j = 1:numel(Events.wincenter)
    tmp = ind.wins{j};
    ind.wins{j} = cat(1, tmp{:});
end

% Sort the cells by position
% sortby = spikes.computeMedianDuringSpikes(Spk.times_spiking, behavior, sortprop, 'quantile_vec', 'lindist');
behavior.blindist = discretize(behavior.lindist, 100);
sortby = spikes.computeOccupancyNormalizedMean(Spk.times_spiking, behavior, 'blindist');

% Initialize the selected struct
selected = struct();
% ---- For Spk -----
numNeurons = numel(Spk.times_spiking);
selected.spike.times_spiking = cell(size(Spk.times_spiking));
for n = progress(1:numNeurons,'Title','Selecting spikes')
    neuronSpikeTimes = Spk.times_spiking{n};
    keepSpikes = false(size(neuronSpikeTimes));
    for i = 1:size(ranges, 2)
        keepSpikes = keepSpikes | (neuronSpikeTimes >= ranges(1, i) & neuronSpikeTimes <= ranges(2, i));
    end
    selected.spike.times_spiking{n} = sort(neuronSpikeTimes(keepSpikes));
    warning off;
    selected.spike.times_spiking{n} = munge.removeDataGaps(selected.spike.times_spiking{n}, ranges, time_gaps, 'gap_thresh', gap_thresh);
    warning on;
end
selected.spike.timeBinMidPoints = munge.removeDataGaps(Spk.timeBinMidPoints(ind.spike), ranges, time_gaps, 'gap_thresh', gap_thresh, 'ploton', false);
selected.spike.spikeCountMatrix = Spk.spikeCountMatrix(:, ind.spike); 
% ---- For efizz -----
selected.efizz.t = munge.removeDataGaps(efizz.t(ind.efizz), ranges, time_gaps, 'gap_thresh', gap_thresh);
fields = {'S1','S2','Cavg','wpli_avg'};
for f = fields
    selected.efizz.(f{1}) = efizz.(f{1})(ind.efizz, :);
end
for f = fieldnames(avg)'
    if use_rescale_avg
        % clamp quantile above below q
        q = 0.01;
        tmp = avg.(f{1})(ind.efizz);
        Q = quantile(tmp, [q, 1-q]); 
        tmp = munge.clamp(tmp, Q(1), Q(2)); 
        selected.efizz.(f{1}) = rescale(tmp, -1, 1);
    else
        selected.efizz.(f{1}) = avg.(f{1})(ind.efizz);
    end
end
disp("efizz" + newline + strjoin(repmat("-", 1, 25)))
disp(selected.efizz)
% --- For Patterns_overall ---
selected.pattern.X_time = ... 
    munge.removeDataGaps(Patterns_overall(end).X_time(ind.pattern), ... 
        ranges, time_gaps, 'gap_thresh', gap_thresh);
selected.pattern.u = Patterns_overall(pattern_overall_ind).cca.u(ind.pattern,:);
selected.pattern.v = Patterns_overall(pattern_overall_ind).cca.v(ind.pattern,:);
selected.pattern.a = Patterns_overall(pattern_overall_ind).cca.a;
selected.pattern.b = Patterns_overall(pattern_overall_ind).cca.b;
disp("Patterns_overall" + newline + strjoin(repmat("-", 1, 25)))
disp(selected.pattern)
selected.behavior_time = munge.removeDataGaps(selected_rows.time, ranges, time_gaps, 'gap_thresh', gap_thresh);
% --- For eeg ---
selected.lfp.time = munge.removeDataGaps(lfp.hpc.time(ind.lfp), ranges, time_gaps, 'gap_thresh', gap_thresh);
for f = lfp_fields
    if use_rescale_avg
        % clamp quantile above below q
        q = 0.01;
        tmp = lfp.hpc.(f{1})(ind.lfp);
        Q = quantile(tmp, [q, 1-q]); 
        tmp = munge.clamp(tmp, Q(1), Q(2)); 
        selected.lfp.hpc.(f{1}) = rescale(tmp, -1, 1);
    else
        selected.lfp.hpc.(f{1}) = lfp.hpc.(f{1}).data(ind.lfp);
    end
end
for f = lfp_fields
    if use_rescale_avg
        % clamp quantile above below q
        q = 0.01;
        tmp = lfp.pfc.(f{1})(ind.lfp);
        Q = quantile(tmp, [q, 1-q]); 
        tmp = munge.clamp(tmp, Q(1), Q(2)); 
        selected.lfp.pfc.(f{1}) = rescale(tmp, -1, 1);
    else
        selected.lfp.pfc.(f{1}) = lfp.pfc.(f{1}).data(ind.lfp);
    end
end
% WARNING: DO ASSERTIONS
% --- For windowed eeg ---
% Create a cell array to store the selected wincenter times for each event type
selected.wins = cell(1, numel(Events.wincenter));
for j = 1:numel(Events.wincenter)
    % Using the indices from ind.wins, extract the corresponding wincenter times
    selected.wins{j} = Events.wincenter{j}(ind.wins{j});
    selected.cellOfWindows{j} = Events.cellOfWindows{j}(ind.wins{j},:);
end
selected.behavior.time = munge.removeDataGaps(selected_rows.time, ranges, time_gaps, 'gap_thresh', gap_thresh);

%% PLOTTING
plots.raw.AllEfizz

%% Behavior plotting
plots.raw.Trajectories

plots.raw.generate_average_set_map(selected, selected_rows, sets_wanna_plot, 'grid_res', 35, 'sgtitlePrepend', animal + " epoch=" + epoch_option + " trajbound=" + trajbound_option);
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps.png"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps.fig"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps.pdf"));
plots.raw.generate_average_set_map(selected, selected_rows, sets_wanna_plot, 'grid_res', 35, 'split_by_reward', true, 'sgtitlePrepend', animal + " epoch=" + epoch_option + " trajbound=" + trajbound_option);
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_splitbyreward.png"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_splitbyreward.fig"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_splitbyreward.pdf"));
plots.raw.generate_average_set_map(selected, selected_rows, sets_wanna_plot, 'grid_res', 35, 'sgtitlePrepend', "STDEV: "  + animal + " epoch=" + epoch_option + " trajbound=" + trajbound_option, 'useRollingStd', true)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_rollingstd.png"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_rollingstd.fig"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_rollingstd.pdf"));
plots.raw.generate_average_set_map(selected, selected_rows, sets_wanna_plot, 'grid_res', 35, 'sgtitlePrepend', "STDEV: " + animal + " epoch=" + epoch_option + " trajbound=" + trajbound_option, 'useRollingStd', true, 'split_by_reward', true);
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_rollingstd_splitbyreward.png"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_rollingstd_splitbyreward.fig"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_averagemaps_rollingstd_splitbyreward.pdf"));

fig(animal + " epoch=" + epoch_option + " trajbound=" + trajbound_option + " corr");clf;
plots.raw.corr
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_corr.png"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_corr.fig"));
saveas(gcf, fullfile(savefolder, animal + "_epoch=" + epoch_option + "_trajbound=" + trajbound_option + "_corr.pdf"));

close all

end % FOR_EPOCH
end % FOR_TRAJBOUND
end % FOR_ANIMAL

!pushover-cli "Finished 2016 plot"
