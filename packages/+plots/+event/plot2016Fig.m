% Tries to display raw data like the 2016 figure
% 
% Inputs:
%   Spk - struct of spike data
%   Events - struct of event data
%   Patterns_overall - struct of CCA data
%   behavior - table of behavior data

%% LOAD DATA

Option = option.defaults();
Option.animal = "ZT2";
Option.preProcess_zscore    = true;
Option.analysis.rankRegress = false;
Option.analysis.cca = false;
Option.save = false;
if ~exist("efizz", "var")
    commsubspaceToPath
    TheScript;
    load(Option.animal + "spectralBehavior.mat");
end
if ~exist('behavior', 'var')
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
end

%% SETTINGS
delta_band = [0.5, 4];
theta_band = [6, 12];
ripple_band = [150, 250];
gap_thresh = 0.5; % Adjust as necessary
use_rescale_avg = true;


%% INDEXING PROPERTIES

% Find frequency indices
delta_idx  = find(efizz.f >= delta_band(1) & efizz.f <= delta_band(2));
theta_idx  = find(efizz.f >= theta_band(1) & efizz.f <= theta_band(2));
ripple_idx = find(efizz.f >= ripple_band(1) & efizz.f <= ripple_band(2));

% Compute average power in the theta and ripple bands for both regions
avg.theta_hpc  = mean(efizz.S1(:, theta_idx), 2);
avg.ripple_hpc = mean(efizz.S1(:, ripple_idx), 2);
avg.theta_pfc  = mean(efizz.S2(:, theta_idx), 2);
avg.ripple_pfc = mean(efizz.S2(:, ripple_idx), 2);

% Define hpc cells
cells.hpc = Spk.areaPerNeuron == "CA1";
cells.pfc = Spk.areaPerNeuron == "PFC";

%% TIME :: Finding Valid Gaps in Data

%% Specify the epoch, traj, trajall, trajclass, trajbound (use [] if not specifying)
epoch_option = 6; % e.g., 'EpochName'
traj_option = []; % a specific trajectory index within an epoch
trajall_option = []; % a specific trajectory index over all epochs
trajclass_option = []; % 1 through 4
trajbound_option = 0; % 0 or 1

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

figure(fig("Time segments"));clf;
tiledlayout('flow'); nexttile
plot(1:size(ranges, 2), ranges(1,:), 'bo', 'DisplayName', 'Start Times');
yline(ranges(1,:), ':', 'Color', 'b');
hold on;
plot(1:size(ranges, 2), ranges(2,:), 'rx', 'DisplayName', 'End Times');
yline(ranges(2,:), ':', 'Color', 'r');
xlabel('Time');
ylabel('Segment Number');
% legend();
indsof = @(x) 1:numel(x);
hold on; plot(Spk.timeBinMidPoints(ind.spike), 'kx', 'DisplayName', 'Spikes');
hold on; plot(selected.spike.timeBinMidPoints, 'gx', 'DisplayName', 'Selected Spikes');
nexttile
plot(Spk.timeBinMidPoints(ind.spike), selected.spike.timeBinMidPoints, 'kx');

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
end
ind.spike = cat(2, ind.spike{:});
ind.efizz = cat(2, ind.efizz{:});
ind.pattern = cat(2, ind.pattern{:});

% Sort the cells by position
spike_pos = spikes.computeMedianDuringSpikes(Spk.times_spiking, behavior, 'lindist');

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
selected.spike.timeBinMidPoints = munge.removeDataGaps(Spk.timeBinMidPoints(ind.spike), ranges, time_gaps, 'gap_thresh', gap_thresh, 'ploton', true);
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
        tmp = avg.(f{1})(ind.efizz)
        Q = quantile(tmp, [q, 1-q]); 
        tmp = munge.clamp(tmp, Q(1), Q(2)); 
        selected.efizz.(f{1}) = rescale(tmp, -1, 1);
    else
        selected.efizz.(f{1}) = avg.(f{1})(ind.efizz)
    end
end
disp("efizz" + newline + strjoin(repmat("-", 1, 25)))
disp(selected.efizz)
% --- For Patterns_overall ---
selected.pattern.X_time = munge.removeDataGaps(Patterns_overall(end).X_time(ind.pattern), ranges, time_gaps, 'gap_thresh', gap_thresh);
selected.pattern.u = Patterns_overall(end).cca.u(ind.pattern,:);
selected.pattern.v = Patterns_overall(end).cca.v(ind.pattern,:);
disp("Patterns_overall" + newline + strjoin(repmat("-", 1, 25)))
disp(selected.pattern)
selected.behavior_time = munge.removeDataGaps(selected_rows.time, ranges, time_gaps, 'gap_thresh', gap_thresh);

%% PLOTTING
fig("2016 Figure");
clf;
width = 600;
height = 1200;
position = [100, 100, width, height];
set(gcf, 'Position', position);

% Define the relative heights
rel_heights = [2, 1, 0.5, 0.5, 0.5, 0.5];
names = ["HPC", "PFC", "Theta", "Ripple", "U", "V"];
normalized_heights = rel_heights / sum(rel_heights);
cumulative_heights = [0, cumsum(normalized_heights)];
% Define user-defined offsets (if needed)
user_offset.theta.hpc = 5; % or whatever value you want
user_offset.theta.pfc = 20; % example value
user_offset.ripple.hpc = -7.5; % or whatever value you want
user_offset.ripple.pfc = 7.5; % example value
% Compute auto offsets based on the middle of the visualization range
visualize.theta = [0, 40];
visualize.ripple = [130, max(efizz.f)];
auto_offset.theta = mean(visualize.theta);
auto_offset.ripple = mean(visualize.ripple);
rowcnt = 0;
rowdict = containers.Map();

ax = gobjects(length(rel_heights), 1); % Pre-allocate memory for axis handles

for i = 1:length(rel_heights)
    % Create a subplot based on calculated positions
    if rel_heights(i) > 0
        rowcnt = rowcnt + 1;
        rowdict(names(i)) = rowcnt;
    else
        continue;
    end
    disp("plotting " + names(i) + " on row " + rowcnt + ", plottype = " + i);
    ax(rowcnt) = subplot('Position', [0.1, cumulative_heights(rowcnt) + 0.01, 0.8, normalized_heights(i) - 0.02]);

    switch i
        case 1
            % Spike Raster for Hippocampal Cells
            plotSpikeRaster(selected.spike.times_spiking(cells.hpc), 'PlotType', 'vertline', 'AutoLabel', false);
            ylabel('HPC Neuron');
            title('Spike Raster of Hippocampal Cells');
        case 2
            % Spike Raster for Prefrontal Cells
            plotSpikeRaster(selected.spike.times_spiking(cells.pfc), 'PlotType', 'vertline', 'AutoLabel', false);
            ylabel('PFC Neuron');
            title('Spike Raster of Prefrontal Cells');
        case 3  % Theta Band
            imagesc(selected.efizz.t, efizz.f, log10(selected.efizz.S1')); set(gca, 'YDir', 'normal', 'YLim', visualize.theta); 
            hold on;
            set(gca, 'Layer', 'top');
            scale = median(structfun(@abs,user_offset.theta));
            plot(selected.efizz.t, scale*selected.efizz.theta_hpc + auto_offset.theta + user_offset.theta.hpc, 'w', 'DisplayName', 'HPC Theta', 'LineWidth', 1.5);
            set(gca, 'Layer', 'top');
            plot(selected.efizz.t, scale*selected.efizz.theta_pfc + auto_offset.theta + user_offset.theta.pfc, 'r', 'DisplayName', 'PFC Theta', 'LineWidth', 1.5);
            ylabel('Frequency (Hz)');
            title('Theta Band Average Power for HPC and PFC');
            legend('Location', 'northwest');
        case 4  % Ripple Band
            tmp = log10(abs(selected.efizz.S1'));
            imagesc(selected.efizz.t, efizz.f, tmp);
            [M,m] = deal(max(tmp(:, ripple_idx),[],'all'), min(tmp(:, ripple_idx),[],'all'));
            clim([m, M]);
            set(gca, 'YDir', 'normal', 'YLim', visualize.ripple);
            hold on;
            scale = 2*median(structfun(@abs,user_offset.ripple));
            plot(selected.efizz.t, scale*selected.efizz.ripple_hpc + auto_offset.ripple + user_offset.ripple.hpc, 'w', 'DisplayName', 'HPC Ripple', 'LineWidth', 1.5);
            plot(selected.efizz.t, scale*selected.efizz.ripple_pfc + auto_offset.ripple + user_offset.ripple.pfc, 'r', 'DisplayName', 'PFC Ripple', 'LineWidth', 1.5);
            ylabel('Frequency (Hz)');
            title('Ripple Band Average Power for HPC and PFC');
            legend('Location', 'northwest');
        case 5  % U components, CCA
            plot(selected.pattern.X_time, selected.pattern.u(:,1), 'DisplayName', 'HPC U');
            hold on;
            plot(selected.pattern.X_time, selected.pattern.u(:,2), 'DisplayName', 'HPC U, 2', 'Color', [0 1 0 0.35], 'LineWidth', 0.5);
            yline(0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
            ylabel('U Component');
            title('U Component of hpc-pfc communication');
        case 6  % V components, CCA
            plot(selected.pattern.X_time, selected.pattern.v(:,1), 'DisplayName', 'PFC V');
            hold on;
            yline(0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
            plot(selected.pattern.X_time, selected.pattern.v(:,2), 'DisplayName', 'PFC V, 2', 'Color', [0 1 0 0.35], 'LineWidth', 0.5);
            ylabel('V Component');
            title('V Component of hpc-pfc communication');
    end
end
linkaxes(ax, 'x');  % Link all the x-axes together

% Overlay Linear Distance on the HPC spike raster
lindist_to_sections = ["HPC", "U", "V"];
for sect = lindist_to_sections
    if rowdict.isKey(sect) % Check if HPC subplot exists
        ax2 = axes('Position', get(ax(rowdict(sect)), 'Position'), 'Color', 'none');
        yyaxis(ax2, 'right');
        hold on;
        % Define the limits for the secondary Y-axis (adjust these as per your needs)
        y2_lim = [min(behavior.lindist), max(behavior.lindist)];
        % ylim(ax2, y2_lim);
        p1 = plot(selected.behavior_time, selected_rows.lindist, 'Color', [1 0 0 0.35], 'Parent', ax2, 'LineWidth', 2);
        % set(p1, 'Color', [1 0 0 0.5]);  % Here, [1 0 0] is red color and 0.35 is the alpha (transparency)
        ylabel('Linear Distance');
        % ax2.YAxisLocation = 'right';  % Ensure the secondary y-axis is on the right
        ax2.Box = 'off';  % This removes the box to make it cleaner
        linkaxes([ax(rowdict("HPC")), ax2], 'x');  % Link the x-axes together
    end
end
set(gcf, 'Position', get(0, 'Screensize'));  % Maximize the figure window
