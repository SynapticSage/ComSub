% Tries to display raw data like the 2016 figure
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
% - window 
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
wins = [1, 3];
shadeOption = false; % show windows : set to true to draw rectangles, and false to skip
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
all_animals = setdiff(all_animals, ["ZT2", "ER1"]);

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
savefolder_opt = "showwin=" + num2str(shadeOption) + "colorbycomp=" + num2str(colorbycomp) + "dim3=" + num2str(dim3) + "sortprop=" + sortprop + "sortdir=" + sortdir;
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
avg.theta_hpc  = mean(efizz.S1(:, theta_idx), 2);
avg.ripple_hpc = mean(efizz.S1(:, ripple_idx), 2);
avg.theta_pfc  = mean(efizz.S2(:, theta_idx), 2);
avg.ripple_pfc = mean(efizz.S2(:, ripple_idx), 2);

% Define hpc cells
cells.hpc = Spk.areaPerNeuron == "CA1";
cells.pfc = Spk.areaPerNeuron == "PFC";

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
selected.lfp.hpc.data = lfp.hpc.data(ind.lfp);
selected.lfp.pfc.data = lfp.pfc.data(ind.lfp);
selected.lfp.hpc.theta = lfp.hpc.theta.data(ind.lfp);
selected.lfp.pfc.theta = lfp.pfc.theta.data(ind.lfp);
selected.lfp.hpc.ripple = lfp.hpc.ripple.data(ind.lfp);
selected.lfp.pfc.ripple = lfp.pfc.ripple.data(ind.lfp);
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
f = fig(animal + " 2016 Figure epoch " + epoch_option + " " + sortprop + " trajbound " + trajbound_option);
figdict = containers.Map();
figdict("figure") = f;
clf;

% Define the relative heights
rel_heights = [2, 1, 0.5, 0.5, 0.5, 0.5];
names = ["HPC", "PFC", "Theta", "Ripple","U", "V"];
normalized_heights = rel_heights / sum(rel_heights);
cumulative_heights = [0, cumsum(normalized_heights)];
% Define user-defined offsets (if needed)
user_offset.theta.hpc  = 5; % or whatever value you want
user_offset.theta.pfc  = 20; % example value
user_offset.ripple.hpc = -7.5; % or whatever value you want
user_offset.ripple.pfc = 7.5; % example value
% Compute auto offsets based on the middle of the visualization range
visualize.theta    = [0, 40];
visualize.ripple   = [130, max(efizz.f)];
auto_offset.theta  = mean(visualize.theta);
auto_offset.ripple = mean(visualize.ripple);

ax = gobjects(length(rel_heights), 1); % Pre-allocate memory for axis handles
rowcnt = 0;
for i = 1:length(rel_heights)
    % Create a subplot based on calculated positions
    if rel_heights(i) > 0
        rowcnt = rowcnt + 1;
        figdict(names(i)) = rowcnt;
    else
        continue;
    end
    disp("plotting " + names(i) + " on row " + rowcnt + ", plottype = " + i);
    figure(figdict("figure"));
    ax(rowcnt) = subplot('Position', [0.1, cumulative_heights(rowcnt) + 0.01, 0.8, normalized_heights(i) - 0.02]);

    switch char(names(i))
        case 'HPC'
            % Spike Raster for Hippocampal Cells
            % figure;tiledlayout('flow');nexttile;imagesc(selected.pattern.a);nexttile;imagesc(sel);
            cla
            spk = selected.spike.times_spiking(cells.hpc);
            [~, isort] = sort(sortby(cells.hpc), sortdir);
            if colorbycomp
                if pattern_overall_ind ~= numel(Patterns_overall)
                    error("colorbycomp only works for pattern_overall_ind = " + numel(Patterns_overall));
                end
                a = selected.pattern.a;
                stds = std(a, [], 1);
                sel = bsxfun(@gt, a, 1*stds) | bsxfun(@lt, a, -1*stds);
                for i = 1:size(spk,2)
                    if isempty(spk{i})
                        continue;
                    end
                    sptmp = spk;
                    [sptmp{1:i-1}] = deal([nan]);
                    [sptmp{i+1:end}] = deal([nan]);
                    color = [0 0; 0 1; 1 0] * sel(i,1:2)';
                    darken = 1.2;
                    color = color ./ (norm(color)*darken);
                    color = fillmissing(color, 'constant', 0);
                    lineformat = struct('Color', color);
                    hold on
                    plotSpikeRaster(sptmp, 'PlotType', 'vertline', 'LineFormat', lineformat, ...
                        'VertSpikePosition', isort(i));
                    hold on;
                end
            else
                plotSpikeRaster(spk(isort), 'PlotType', 'vertline', 'AutoLabel', false);
                get(gca,'Children');
            end
            ylabel('HPC Neuron');
            title('Spike Raster of Hippocampal Cells');
        case 'PFC'
            % Spike Raster for Prefrontal Cells
            spk = selected.spike.times_spiking(cells.pfc);
            [~, isort] = sort(sortby(cells.pfc), sortdir);
            if colorbycomp
                if pattern_overall_ind ~= numel(Patterns_overall)
                    error("colorbycomp only works for pattern_overall_ind = " + numel(Patterns_overall));
                end
                a = selected.pattern.a;
                stds = std(a, [], 1);
                sel = bsxfun(@gt, a, 1*stds) | bsxfun(@lt, a, -1*stds);
                for i = 1:size(spk,2)
                    if isempty(spk{i})
                        continue;
                    end
                    sptmp = spk;
                    [sptmp{1:i-1}] = deal([nan]);
                    [sptmp{i+1:end}] = deal([nan]);
                    color = [0 0; 0 1; 1 0] * sel(i,1:2)';
                    darken = 1.2;
                    color = color ./ (norm(color)*darken);
                    color = fillmissing(color, 'constant', 0);
                    lineformat = struct('Color', color);
                    hold on
                    plotSpikeRaster(sptmp, 'PlotType', 'vertline', 'LineFormat', lineformat, ...
                        'VertSpikePosition', isort(i));
                    hold on;
                end
            else
                plotSpikeRaster(spk(isort), 'PlotType', 'vertline', 'AutoLabel', false);
                get(gca,'Children');
            end
            ylabel('PFC Neuron');
            title('Spike Raster of Prefrontal Cells');
        case 'Theta'
            tmp = log10(abs(selected.efizz.S1'));
            imagesc(selected.efizz.t, efizz.f, tmp); 
            set(gca, 'YDir', 'normal', 'YLim', visualize.theta); 
            [M,m] = deal(quantile(tmp(:, theta_idx),0.95,'all'), quantile(tmp(:, theta_idx),0.05,'all'));
            cmocean('thermal');
            clim([m, M]);
            hold on;
            set(gca, 'Layer', 'top');
            scale = median(structfun(@abs,user_offset.theta));
            plot(selected.efizz.t, scale*selected.efizz.theta_hpc + auto_offset.theta + user_offset.theta.hpc, 'w', 'DisplayName', 'HPC Theta', 'LineWidth', 1.5);
            set(gca, 'Layer', 'top');
            plot(selected.efizz.t, scale*selected.efizz.theta_pfc + auto_offset.theta + user_offset.theta.pfc, 'r', 'DisplayName', 'PFC Theta', 'LineWidth', 1.5);
            ylabel('Frequency (Hz)');
            title('Theta Band Average Power for HPC and PFC');
            legend('Location', 'northwest');
        case 'Ripple'
            tmp = log10(abs(selected.efizz.S1'));
            imagesc(selected.efizz.t, efizz.f, tmp);
            [M,m] = deal(quantile(tmp(:, ripple_idx),0.7,'all'), quantile(tmp(:, ripple_idx),0.01,'all'));
            cmocean('thermal');
            clim([m, M]);
            set(gca, 'YDir', 'normal', 'YLim', visualize.ripple);
            hold on;
            scale = 2*median(structfun(@abs,user_offset.ripple));
            plot(selected.efizz.t, scale*selected.efizz.ripple_hpc + auto_offset.ripple + user_offset.ripple.hpc, 'w', 'DisplayName', 'HPC Ripple', 'LineWidth', 1.5);
            plot(selected.efizz.t, scale*selected.efizz.ripple_pfc + auto_offset.ripple + user_offset.ripple.pfc, 'r', 'DisplayName', 'PFC Ripple', 'LineWidth', 1.5);
            ylabel('Frequency (Hz)');
            title('Ripple Band Average Power for HPC and PFC');
        case 'U'
            % Shading for U Component 1
            % fill(X, Y1, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            hold on;
            plot(selected.pattern.X_time, selected.pattern.u(:,1), 'DisplayName', 'HPC U');
            plots.fill_curve(selected.pattern.X_time, selected.pattern.u(:,1), 'b')
            alpha(0.5);

            % Shading for U Component 2
            % fill(X, Y2, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(selected.pattern.X_time, selected.pattern.u(:,2), 'DisplayName', 'HPC U, 2', 'Color', [0 1 0 0.35], 'LineWidth', 0.5);
            plots.fill_curve(selected.pattern.X_time, selected.pattern.u(:,2), 'g')
            alpha(0.3);
            if dim3
                % Shading for U Component 3
                % fill(X, Y3, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                plot(selected.pattern.X_time, selected.pattern.u(:,dim3), 'DisplayName', 'HPC U, 3', 'Color', [1 0 0 0.2], 'LineWidth', 0.5);
                plots.fill_curve(selected.pattern.X_time, selected.pattern.u(:,3), 'r')
                alpha(0.1);
            end
            yline(0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
            ylabel('U Component');
            title('U Component of hpc-pfc communication');

        case 'V'
            % Shading for V Component 1
            plot(selected.pattern.X_time, selected.pattern.v(:,1), 'DisplayName', 'PFC V');
            hold on;
            plots.fill_curve(selected.pattern.X_time, selected.pattern.v(:,1), 'b');
            alpha(0.5);
            if dim3
                % Shading for V Component 3
                plot(selected.pattern.X_time, selected.pattern.v(:,dim3), 'DisplayName', 'PFC V, 3', 'Color', [1 0 0 0.2], 'LineWidth', 0.5);
                plots.fill_curve(selected.pattern.X_time, selected.pattern.v(:,3), 'r');
                alpha(0.1);
            end
            % Shading for V Component 2
            yline(0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
            plot(selected.pattern.X_time, selected.pattern.v(:,2), 'DisplayName', 'PFC V, 2', 'Color', [0 1 0 0.35], 'LineWidth', 0.5);
            plots.fill_curve(selected.pattern.X_time, selected.pattern.v(:,2), 'g');
            alpha(0.3);
            ylabel('V Component');
            title('V Component of hpc-pfc communication');

        case 'Raw'

            plot(lfp.hpc.time, lfp.hpc.data(:,1));
            hold on
            plot(lfp.pfc.time, lfp.pfc.data(:,1));
            title('Raw Voltage');
            ylabel('Amplitude');
        
        case 'RawTheta'

            plot(lfp.hpc.theta.time, lfp.hpc.theta.data(:,1));
            hold on
            plot(lfp.pfc.theta.time, lfp.pfc.theta.data(:,1));
            title('Theta Power');
            ylabel('Amplitude');

        case 'RawRipple'

            plot(lfp.hpc.theta.time, lfp.hpc.theta.data(:,1));
            plot(lfp.pfc.theta.time, lfp.pfc.theta.data(:,1));
            title('Ripple Power');
            ylabel('Amplitude');

    end

    % Network pattern
    if ~isempty(shadeOption)
        % Get the current axis
        thisax = gca;
       % Make the current axis active
        hold on;
        % Iterate over the windows in cellOfWindows
        for netpat = 1:numel(selected.cellOfWindows)
        if ismember(netpat, shadeOption)
        for win_idx = 1:size(selected.cellOfWindows{netpat}, 1)
            currWin = selected.cellOfWindows{netpat}(win_idx, :);
            % Shade the region using patch
            yLimits = ylim(thisax);
            patch('XData', [currWin(1), currWin(1), currWin(2), currWin(2)], ...
                  'YData', [yLimits(1), yLimits(2), yLimits(2), yLimits(1)], ...
                  'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none', ...
                  'FaceAlpha', 0.2, 'Parent', ax(rowcnt));
            hold on;
        end
        end
        end
        % Restore the hold state of the axis
        hold off; 
    end

end
linkaxes(ax, 'x');  % Link all the x-axes together

% Overlay Linear Distance on the HPC spike raster
lindist_to_sections = ["HPC", "U", "V"];
for sect = lindist_to_sections
    if figdict.isKey(sect) % Check if HPC subplot exists
        % figure(figdict(sect));
        ax2 = axes('Position', get(ax(figdict(sect)), 'Position'), 'Color', 'none');
        yyaxis(ax2, 'right');
        hold on;
        % Define the limits for the secondary Y-axis (adjust these as per your needs)
        y2_lim = [min(behavior.lindist), max(behavior.lindist)];
        ylim(ax2, y2_lim);
        p1 = plot(selected.behavior_time, selected_rows.lindist,...
        'Color', [0 0 0 0.35], 'Parent', ax2, 'LineWidth', 3.5);
        % set(p1, 'Color', [1 0 0 0.5]);  % Here, [1 0 0] is red color and 0.35 is the alpha (transparency)
        ylabel('Linear Distance');
        % ax2.YAxisLocation = 'right';  % Ensure the secondary y-axis is on the right
        ax2.Box = 'off';  % This removes the box to make it cleaner
        linkaxes([ax(figdict("HPC")), ax2], 'x');  % Link the x-axes together
    end
end
sgt = gcf;
set(gcf, 'Position', get(0, 'Screensize'));  % Maximize the figure window
r = @(x) string(replace(replace(replace(x, " ", "_"), ":", "_"), newline, "_"));
saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.png'));
saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.fig'));

% width = 600;
% height = 1200;
% position = [100, 100, width, height];
% set(gcf, 'Position', position);

%% Behavior plotting

% Extract unique trajectories
unique_traj = unique(selected_rows.trajall);
% Randomly sample 10,000 points from the main behavior table
rand_indices = randperm(height(behavior), min(10000, height(behavior)));
sampled_X = behavior.X(rand_indices);
sampled_Y = behavior.Y(rand_indices);
xlims = quantile(sampled_X, [0.01, 0.99]);
ylims = quantile(sampled_Y, [0.01, 0.99]);

% Parameters for the behavior plot
for set_to_plot = sets_wanna_plot(:)' %FOR_SET_TO_PLOT

    % time_name = "X_time";
    % data_name = "v";
    % data_col  = 2;
    time_name = set_to_plot{1}(3);
    data_name = set_to_plot{1}(1);
    data_col  = double(set_to_plot{1}(2));
    normalized = "minmaxabs";
    pick = 12;
    L = length(unique_traj);
    randset = sort(randperm(L, min(L,pick)));
    background_color = [0 0 0];
    % background_color = [1 1 1];

    % Create figure and tiled layout
    f = fig(animal + "Behavior Trajectories " +  + epoch_option + " " + sortprop + " trajbound " + trajbound_option + newline + " dataname: " + data_name + ", col: " + data_col);clf;
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

    % Color data based on some field in the selected structure
    data_for_coloring = interp1(selected.pattern.(time_name), selected.pattern.(data_name)(:,data_col), selected.behavior.time, 'nearest','extrap');
    extrap_locations = isnan(interp1(selected.pattern.(time_name), selected.pattern.(data_name)(:,data_col), selected.behavior.time, 'nearest'));
    % Colormap for the data (can be changed to any other colormap)
    cm=cmocean('balance', 1024);
    % Normalize the interpolated data to be within [0, 1]
    if normalized == "minmax"
        normalized_data = (data_for_coloring - min(data_for_coloring)) / (max(data_for_coloring) - min(data_for_coloring));
        clims = [min(data_for_coloring), max(data_for_coloring)];
    elseif normalized == "minmaxabs"
        normalized_data = ((data_for_coloring) / (max(abs(data_for_coloring))))/2 + 0.5
        clims = [-max(abs(data_for_coloring)), max(abs(data_for_coloring))];
    end
    % Map the normalized data to colormap indices
    color_indices = round(normalized_data * (size(cm, 1) - 1) + 1);
    % Convert indices to RGB values
    colors_for_plotting = cm(color_indices, :);
    colors_for_plotting(extrap_locations, :) = repmat([0 0 0], sum(extrap_locations), 1);

    % Loop through each unique trajectory and plot X and Y positions
    for traj_num = progress(unique_traj(randset)','Title', 'Plotting Trajectories')
        % Extract the current trajectory data
        inds_traj = selected_rows.trajall == traj_num;
        curr_traj = selected_rows(inds_traj, :);
        if ~isempty(data_for_coloring) && any(extrap_locations(inds_traj))
            % If there are any extrapolated locations, plot them in black
            warning('There are extrapolated locations in trajectory %d', traj_num);
            continue;
        end
        
        % Next tile for the current trajectory
        ax = nexttile;
        
        % Plot the sampled background trajectory
        plot(ax, sampled_X, sampled_Y, '.', 'Color', [0.5 0.5 0.5]);
        hold on;
        
        if ~isempty(data_for_coloring)
            % Plot the trajectory colored by the data
            scatter(ax, curr_traj.X, curr_traj.Y, 40, colors_for_plotting(inds_traj,:), 'filled');
            colormap(ax, cm);
            h=colorbar;
            h.Label.String = 'u';
            h.TickLabels = cellstr(string(linspace(clims(1), clims(2), 5)'));
            h.Ticks = linspace(0, 1, 5);
        else
            % Plot the trajectory without coloring
            plot(ax, curr_traj.X, curr_traj.Y, 'k');
        end
        
        % Label for the start and end of trajectory
        plot(ax, curr_traj.X(1), curr_traj.Y(1), 'go');
        plot(ax, curr_traj.X(end), curr_traj.Y(end), 'ro');
        
        % Set title
        if curr_traj.rewarded(1) == 1
            rew_str = 'Rewarded';
        else
            rew_str = 'Not Rewarded';
        end
        if curr_traj.leftright(1)
            lr_str = 'Left';
        else
            lr_str = 'Right';
        end
        title(['Epoch:' num2str(epoch_option) ' Traj: ' num2str(traj_num) newline 'TrajBound: ', num2str(curr_traj.trajbound(1)), ', ', rew_str, ', ', lr_str]);
        
        % Additional axis properties (if needed)
        xlabel('X Position');
        ylabel('Y Position');
        set(ax, 'Color', background_color);
        axis equal; % Make the X and Y axis scales the same
        ylim(ylims);
        xlim(xlims);
        grid on;
    end
    sgtitle("Behavior Trajectories " +  + epoch_option + " " + sortprop + " trajbound " + trajbound_option + newline + " dataname: " + data_name + ", col: " + data_col);
    sgt = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));  % Maximize the figure window
    r = @(x) string(replace(replace(replace(x, " ", "_"), ":", "_"),newline,"_"));
    saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.png'));
    saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.pdf'));
    saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.fig'));

end % FOR_EPOCH
end % FOR_TRAJBOUND
end % FOR_SET_TO_PLOT
end % FOR_ANIMAL

!pushover-cli "Finished 2016 plot"
