function out = triggered_spectrogram(Patterns_overall, Spk, efizz, varargin)
% SPEC_ANALYSIS Get windows of efizz data around comm sub threshold crossings
%   
% Inputs:
%   Patterns_overall: struct array containing the CCA patterns for each
%       comm sub
%   Spk: struct containing the spike data
%   efizz: struct containing the efizz data
%
% Outputs:
%   out: struct containing the average efizz data for each comm sub

ip = inputParser();
ip.addParameter('ploton', false, @islogical);
ip.addParameter('windowsize', 50, @isnumeric);
ip.addParameter('specNames', {'S1','S2','Cavg','wpli_avg','phi'}, @iscellstr);
ip.addParameter('freq_ylims', [0, 50], @isnumeric);
ip.addParameter('components', [1,2,3,4,5], @isnumeric);
ip.addParameter('figAppend', "");
ip.addParameter('folder', "triggered_spectrogram", @(x) isstring(x));
ip.addParameter('quantile_threshold', 0.95, @isnumeric);
ip.addParameter('runtype', 1, @isnumeric); % 1 = run, 0 = rest
ip.addParameter('boots', 0, @isnumeric);
ip.addParameter('based_on', "mean");
ip.parse(varargin{:});
Opt = ip.Results;
Opt.folder = string(Opt.folder) + "_run=" + Opt.runtype + filesep;
boots = Opt.boots;

if isstruct(Opt.figAppend)
    scalar_info = Opt.figAppend;
    tmp = string(struct2cell(scalar_info));
    tmp = strjoin(tmp, "_");
    Opt.figAppend = tmp;
else
    scalar_info = struct();
end

if ~exist(figuredefine(Opt.folder), 'dir')
    mkdir(figuredefine(Opt.folder));
end
if ~isempty(Opt.folder) && ~startsWith(Opt.folder, filesep)
    Opt.figAppend = "_" + Opt.figAppend;
end
% Choose your quantile threshold
quantile_threshold = Opt.quantile_threshold;
% Define the window size around the threshold crossing (e.g., 10 time bins
% before and after)
window_size = Opt.windowsize;
% Determine total window size
total_window_size = 2 * window_size + 1;  % accounting for the point itself
% Find indices of pattern with name="overall"
% overall_ind = find(isequal([Patterns_overall.name], "Overall"));

disp("Starting with " + length(Patterns_overall) + " patterns and " + ...
    length(Opt.components) + " components and ploton = " + Opt.ploton);

Opt.means = [];
Opt.mins = [];
for field = Opt.specNames
    if contains(field, "wpli")
        efizz.(field{:})(isnan(efizz.(field{:}))) = 0;
    end
    Opt.means.(field{:}) = mean(efizz.(field{:}), 1);
    Opt.mins.(field{:}) = min(efizz.(field{:}), [], 1);
end

% Loop over all patterns
for i = progress(1:numel(Patterns_overall), 'Title', 'Patterns')

    [d, p] = ind2sub(size(Patterns_overall), i);

    if isempty(Patterns_overall(i).cca) || isempty(Patterns_overall(i).cca.u) || ...
       isempty(Patterns_overall(i).cca.v) || ...
        Patterns_overall(i).directionality == "hpc-hpc" 
        continue
    end

    % Extract u and v from the current pattern
    runs = Spk.sessionTypePerBin==Opt.runtype;
    a = Patterns_overall(i).cca.a;
    b = Patterns_overall(i).cca.b;
    ca1 = Spk.areaPerNeuron == "CA1";
    pfc = Spk.areaPerNeuron == "PFC";
    source = Spk.spikeRateMatrix(ca1, runs);
    target = Spk.spikeRateMatrix(pfc, runs);
    epochs = Spk.epochPerBin(runs);
    u = source' * a;
    v = target' * b;
    sp_time = Spk.timeBinMidPoints(runs);
    assert(numel(sp_time) == size(u,1), ... 
        'Spike time and u must be the same length');

    ustd = movstd(u, 7);
    vstd = movstd(v, 7);

    % Calculate the quantile threshold for u and v
    u_threshold = quantile(u, quantile_threshold);
    v_threshold = quantile(v, quantile_threshold);
    u_threshold_lower = quantile(u, 1-quantile_threshold);
    v_threshold_lower = quantile(v, 1-quantile_threshold);

    % Find the time bins where u or v cross their respective thresholds
    if Opt.based_on == "mean"
        u_above_threshold = u > u_threshold;
        v_above_threshold = v > v_threshold;
        all_threshold_crossed = u_above_threshold & v_above_threshold;
    elseif Opt.based_on == "mean_u_above"
        u_above_threshold = u > u_threshold;
        all_threshold_crossed = u_above_threshold;
    elseif Opt.based_on == "mean_v_above"
        v_above_threshold = v > v_threshold;
        all_threshold_crossed = v_above_threshold;
    elseif Opt.based_on == "mean_u_below"
        u_below_threshold = u < u_threshold_lower;
        all_threshold_crossed = u_below_threshold;
    elseif Opt.based_on == "mean_v_below"
        v_below_threshold = v < v_threshold_lower;
        all_threshold_crossed = v_below_threshold;
    elseif Opt.based_on == "std"
        u_stdthresh           = ustd > quantile(ustd, quantile_threshold);
        v_stdthersh           = vstd > quantile(vstd, quantile_threshold);
        all_threshold_crossed = u_stdthresh & v_stdthersh;
    elseif Opt.based_on == "std_u"
        u_stdthresh           = ustd > quantile(ustd, quantile_threshold);
        all_threshold_crossed = u_stdthresh;
    elseif Opt.based_on == "std_v"
        v_stdthersh           = vstd > quantile(vstd, quantile_threshold);
        all_threshold_crossed = v_stdthersh;
    elseif Opt.based_on == "mean_and_std"
        u_above_threshold = u > u_threshold;
        u_stdthresh       = ustd > quantile(ustd, quantile_threshold);
        v_above_threshold = v > v_threshold;
        v_stdthersh       = vstd > quantile(vstd, quantile_threshold);
        all_threshold_crossed = (u_above_threshold & v_above_threshold) & ...
                                (u_stdthresh & v_stdthersh);
    else
        error("Invalid value for Opt.based_on: " + Opt.based_on);
    end

    if isempty(Patterns_overall(i).nameFull)
        name = Patterns_overall(i).name;
    else
        name = Patterns_overall(i).nameFull;
    end

    for comp = progress(Opt.components(:)','Title','Components')

        threshold_crossed = all_threshold_crossed(:,comp);  % make sure it's a column vector
        if ~any(threshold_crossed)
            continue
        end

        % Determine the corresponding times in efizz
        threshold_crossed_times = sp_time(threshold_crossed); 

        % Find the indices of these times in the efizz data
        efizz_indices = interp1(efizz.t, 1:length(efizz.t), threshold_crossed_times, 'nearest');
        spike_indices = find(threshold_crossed);
        efizz_indices = efizz_indices(:);
        spike_indices = spike_indices(:);
        both = [efizz_indices, spike_indices];
        both = both(~isnan(both(:,1)),:);
        efizz_indices = both(:,1);
        % spike_indices = both(:,2);
        % Find start/stop times for spike windows
        inds_starts = efizz_indices-window_size;
        inds_stops  = efizz_indices+window_size;
        % Throw away out of bounds indices
        both = both(inds_starts > 0 & inds_stops <= length(efizz.t),:);
        efizz_indices = both(:,1);
        spike_indices = both(:,2);
        % Find start/stop times for spike windows
        sp_start_indices = interp1(sp_time, 1:numel(sp_time), efizz.t(efizz_indices+window_size), 'nearest');
        sp_start_indices = sp_start_indices(:);
        sp_stop_indices  = interp1(sp_time, 1:numel(sp_time), efizz.t(efizz_indices+window_size), 'nearest');
        sp_stop_indices  = sp_stop_indices(:);
        if isempty(sp_start_indices) || isempty(sp_stop_indices)
            continue
        end
        sp_win_size      = mode(abs([sp_stop_indices - spike_indices; spike_indices - sp_start_indices]))+1;
        total_sp_win_size = 2 * sp_win_size + 1;

        % Initialize empty matrices to store the data segments for each efizz
        % variable
        freqsize = size(efizz.S1, 2);
        segments = cell(1, length(Opt.specNames));
        for j = 1:length(Opt.specNames)
            segments{j} = nan(total_window_size, freqsize, length(efizz_indices));
        end
        time_segments    = nan(total_window_size, length(efizz_indices));
        u_segments       = nan(total_sp_win_size, size(u, 2), length(efizz_indices));
        v_segments       = nan(total_sp_win_size, size(v, 2), length(efizz_indices));
        sp_time_segments = nan(total_sp_win_size, length(efizz_indices));
        epoch_segments   = nan(total_sp_win_size, length(efizz_indices));

        % For each efizz index, grab the data in a window around that index
        % for j = progress(1:min(size(both,1),10),'Title','Extracting efizz data')
        for j = progress(1:size(both,1),'Title','Extracting efizz data')
            index = efizz_indices(j);
            spike_index = spike_indices(j);
            eftimes = efizz.t(index-window_size:index+window_size);
            if any(diff(eftimes) > 2)
                continue
            end
            if index-window_size > 0 && index+window_size <= size(efizz.S1, 1)  ...
                && spike_index-sp_win_size > 0 && spike_index+sp_win_size <= size(u, 1)
                % Make sure the window does not exceed the matrix dimensions
                for k = 1:length(Opt.specNames)
                    segments{k}(:,:,j) = ...
                        efizz.(Opt.specNames{k})(index-window_size:index+window_size, :);
                end
                time_segments(:,j) = linspace(efizz.t(index-window_size), ...
                                              efizz.t(index+window_size), ...
                                              total_window_size);
                % Spiking
                sp_time_segments(:,j) = ...
                    sp_time(spike_index-sp_win_size:spike_index+sp_win_size);
                u_segments(:,:,j) = ...
                    u(spike_index-sp_win_size:spike_index+sp_win_size, :);
                v_segments(:,:,j) = ...
                    v(spike_index-sp_win_size:spike_index+sp_win_size, :);
                epoch_segments(:,j) = ...
                    epochs(spike_index-sp_win_size:spike_index+sp_win_size);
            end
        end

        disp("Averaging " + length(efizz_indices) + " segments for " + name + ...
             " direction " + d + " comp " + comp);
        for j = 1:length(Opt.specNames)
            s=segments{j};
            s=s(:,:,~isnan(s(1,1,:)));
            disp("...stats for " + Opt.specNames{j});
            disp("")
            if Opt.specNames{j} == "phi"
                spec_avg.(Opt.specNames{j})      = angle(mean(exp(1i*s), 3));
                spec_stderr.(Opt.specNames{j})   = std(exp(1i*s), 0, 3) / sqrt(size(s, 3));
                if boots
                    disp("...calculating boots=" + boots + " confidence intervals for " + Opt.specNames{j});
                    spec_ci_upper.(Opt.specNames{j}) = bootci(boots, @(x) angle(mean(exp(1i*x), 1)), s);
                    spec_ci_lower.(Opt.specNames{j}) = bootci(boots, @(x) angle(mean(exp(1i*x), 1)), s);
                end
            else
                spec_avg.(Opt.specNames{j})      = mean(s, 3);
                spec_stderr.(Opt.specNames{j})   = std(s, 0, 3) / sqrt(size(s, 3));
                if boots
                    disp("...calculating boots=" + boots + " confidence intervals for " + Opt.specNames{j});
                    spec_ci_upper.(Opt.specNames{j}) = bootci(boots, @(x) mean(x, 1), s);
                    spec_ci_lower.(Opt.specNames{j}) = bootci(boots, @(x) mean(x, 1), s);
                end
            end
        end
        time_segments = time_segments(:,~isnan(time_segments(1,:)));
        time_avg      = mean(time_segments, 2);

        % Length
        u_segments = u_segments(:,:,~isnan(u_segments(1,1,:)));
        v_segments = v_segments(:,:,~isnan(v_segments(1,1,:)));
        % Interpolate UV to match efizz window segment time size
        told = linspace(time_avg(1), time_avg(end), size(u_segments, 1));
        tnew = linspace(time_avg(1), time_avg(end), total_window_size);
        u_average  = mean(u_segments, 3);
        u_average  = interp1(told(:), u_average, tnew, 'linear');
        v_average  = mean(v_segments, 3);
        v_average  = interp1(told(:), v_average, tnew, 'linear');
        u_stderr   = std(u_segments, 0, 3) / sqrt(size(u_segments, 3));
        u_stderr   = interp1(told(:), u_stderr, tnew, 'linear');
        v_stderr   = std(v_segments, 0, 3) / sqrt(size(v_segments, 3));
        v_stderr   = interp1(told(:), v_stderr, tnew, 'linear');
        if boots
            disp("bootstrapping boot=" + boots + " ci...")
            u_ci_upper = bootci(boots, @(x) mean(x, 1), u_segments);
            v_ci_upper = bootci(boots, @(x) mean(x, 1), v_segments);
            disp("...done!")
        end

        out(i,comp).spec_avg = spec_avg;
        out(i,comp).spec_stderr = spec_stderr;
        if boots
            out(i,comp).spec_ci_upper = spec_ci_upper;
            out(i,comp).spec_ci_lower = spec_ci_lower;
        end
        out(i,comp).threshold_crossed_times = threshold_crossed_times;
        if boots
            out(i,comp).u_ci_upper = u_ci_upper;
            out(i,comp).v_ci_upper = v_ci_upper;
        end
        out(i,comp).name = name;
        out(i,comp).comp = comp;
        out(i,comp).direction = d;
        out(i,comp).u_average = u_average;
        out(i,comp).v_average = v_average;
        out(i,comp).u_stderr = u_stderr;
        out(i,comp).v_stderr = v_stderr;
        out(i,comp).u_threshold = u_threshold(comp);
        out(i,comp).v_threshold = v_threshold(comp);
        out(i,comp).epoch = mode(epoch_segments, 1);
        out(i,comp).time_avg = tnew(:);

    end

    disp("Writing parquet files for " + name + " direction " + d);

    % save to table
    try
    [t_uv, t_spec] = table.analyses.trigspec(out, [], efizz.f, scalar_info);
    tablefolder = figuredefine("tables");
    parquetwrite(fullfile(tablefolder, "triggeredspec_uv"+Opt.figAppend+".parquet"), t_uv);
    parquetwrite(fullfile(tablefolder, "triggeredspec_spec"+Opt.figAppend+".parquet"), t_spec);
    catch
        disp("Error writing parquet files for " + name + " direction " + d);
    end

    if Opt.ploton
        try
        disp("Plotting " + name + " direction " + d);
        plots.triggered_spec_struct(out(i,:), efizz, Opt);
        catch
            disp("Error plotting " + name + " direction " + d);
        end
    end
end
