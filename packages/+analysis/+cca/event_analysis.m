function out = event_analysis(Patterns_overall, Spk, Events, Option, behavior, varargin)
%EVENT_ANALYSIS   Calculate CCA r-values for each event in a given pattern.
%
%  out = event_analysis(Components, Patterns_overall, Spk, Events)
%
%  INPUTS:
%  Patterns_overall
%  Spk:  Spike data structure 
%  Events:  Structure containing event times
%  Option:  Structure containing options
%   varargin:  'precompute', 0 or 1 (default 0)
%              'method', 'zscore' or 'spikerate' (default 'zscore')
%
%  OUTPUTS:
%  out:  Structure containing CCA r-values for each event in each pattern

ip = inputParser();
ip.addParameter('precompute', 0, @isscalar);
ip.addParameter('method', 'zscore', @ischar);
ip.addParameter('scalar_struct', []);
ip.addParameter('N', 5, @isscalar);
ip.parse(varargin{:});
Opt = ip.Results;

if isempty(Opt.scalar_struct)
    Opt.scalar_struct = struct(...
        'animal', Option.animal,...
        'genH', Option.genH_name,...
        'zscore', Option.preProcess_zscore);
end

% Assuming 'area1' and 'area2' are the indices of the areas you're interested in
area1 = find(strcmp(Spk.areaPerNeuron, 'CA1'));
area2 = find(strcmp(Spk.areaPerNeuron, 'PFC'));

% Create an empty matrix to store the CCA r-values for each event
out = struct('W', []);
% out = repmat(out, size(Patterns_overall));

szPatterns = size(Patterns_overall);

% Loop over all patterns
for i = progress(1:numel(Patterns_overall), 'Title', 'Event analysis')

    event_u = cell(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_v = cell(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_r_values = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), Opt.N);
    event_u_values = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), Opt.N);
    event_v_values = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), Opt.N);
    event_u_values_mean = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), Opt.N);
    event_v_values_mean = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), Opt.N);
    event_time = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_epoch = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_lindist = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_vel = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_trajbound = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_correct = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    % uv_components = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));

    if isempty(Patterns_overall(i).cca)
        continue;
    end

    % Extract a and b from the current pattern
    if isfield(Patterns_overall(i).cca, 'a') && Opt.precompute >= 1
        disp("Using precomputed...")
        a = Patterns_overall(i).cca.a;
        b = Patterns_overall(i).cca.b;
    else % if not precomputed
        disp("Computing CCA")
        source = Patterns_overall(i).X_source;
        target = Patterns_overall(i).X_target;
        index_source = Patterns_overall(i).index_source;
        index_target = Patterns_overall(i).index_target;
        if Opt.method == "zscore" && ~munge.detectZscore(source)
            % zscore the fioring
            fprintf("...zscore\n")
            source = zscore(source, 0, 2);
            target = zscore(target, 0, 2);
        elseif Opt.method == "spikerate" && munge.detectZscore(source)
            fprintf("...un-zscore\n")
            % undo z-score normalization             
            muFR   = Spk.muFR(Spk.hpc.to_original(index_source));
            stdFR  = Spk.stdFR(Spk.hpc.to_original(index_source));
            source = source .* stdFR + muFR;
            muFR   = Spk.muFR(Spk.pfc.to_original(index_target));
            stdFR  = Spk.stdFR(Spk.pfc.to_original(index_target));
            target = target .* stdFR + muFR;
        end
        [a,b] = canoncorr(source', target');
    end

    [i1,i2] = ind2sub(szPatterns, i);
    i = {i1, i2};
    if i2 <= Option.nPatternAndControl
        p = i2;
    else
        p = 1:Option.nPatternAndControl;
    end
    directionality = Patterns_overall(i{:}).directionality;

    % Loop over all events
    ecnt=0;
    for w = p % uv components
        windows = Events.cellOfWindows{w};
        disp("length windows = " + length(windows))
        for j = 1:length(windows) % events

            % Find the time bins that correspond to the current event
            event_time_bins = find(Spk.timeBinMidPoints >= windows(j,1) &...
                Spk.timeBinMidPoints <= windows(j,2));
            if isempty(event_time_bins)
                ecnt=ecnt+1;
                continue;
            end
            center_time = mean(windows(j,:));

            % Extract the spike data for the two areas during this event
            if directionality == "hpc-pfc"
                area1_spikes = Spk.spikeCountMatrix(area1, event_time_bins);
                area2_spikes = Spk.spikeCountMatrix(area2, event_time_bins);
            elseif directionality == "hpc-hpc"
                continue;
            end

            % Project the spike data onto the space defined by a and b to get u and v
            u = area1_spikes' * a(:,1:Opt.N);
            v = area2_spikes' * b(:,1:Opt.N);

            % Calculate the correlation between u and v
            r = nan(Opt.N,1);
            for n = 1:Opt.N
                r(n) = corr(u(:,n), v(:,n));
            end

            behtimes = interp1(behavior.time, 1:numel(behavior.time), center_time, 'nearest');
            if all(~isnan(behtimes))
                epoch     = behavior.epoch(behtimes);
                lindist   = behavior.lindist(behtimes);
                vel       = behavior.vel(behtimes);
                trajbound = behavior.trajbound(behtimes);
                correct   = behavior.rewarded(behtimes);
            else
                epoch     = nan;
                lindist   = nan;
                vel       = nan;
                trajbound = nan;
                correct   = nan;
            end
            
            % Store the CCA r-value and canonical variates for this event
            event_r_values(w,j,:) = r;  % assuming we are interested in the first canonical correlation
            event_u_values(w,j,:) = mean(u,1);
            event_v_values(w,j,:) = mean(v,1);
            event_u_values_mean(w,j,:)  = repmat(mean(mean(u,1)),1,Opt.N);
            event_v_values_mean(w,j,:)  = repmat(mean(mean(v,1)),1,Opt.N);
            % figure out how many nans
            % arrayfun(@(x) sum(isnan(event_v_values(:,x,:)),'all'), 1:size(event_v_values,2))
            event_u{w,j} = u;
            event_v{w,j} = v;
            event_time(w,j) = center_time;
            event_epoch(w,j) = epoch;
            event_lindist(w,j) = lindist;
            event_vel(w,j) = vel;
            event_trajbound(w,j) = trajbound;
            event_correct(w,j) = correct;
        end
    end
    disp("...done")
    disp("fraction of empty events: " + num2str(ecnt/numel(Events.cellOfWindows)))
    % Store the CCA r-values and canonical variates for this pattern
    out(i{:}).event_r_values = event_r_values;
    out(i{:}).event_u_values = event_u_values;
    out(i{:}).event_v_values = event_v_values;
    out(i{:}).event_u        = event_u;
    out(i{:}).event_v        = event_v;
    out(i{:}).event_time     = event_time;
    out(i{:}).epoch          = event_epoch;
    out(i{:}).lindist        = event_lindist;
    out(i{:}).vel            = event_vel;
    out(i{:}).trajbound      = event_trajbound;
    out(i{:}).correct        = event_correct;
    % out(i{:}).uv_comp        = uv_components;
end

tabappend = struct2cell(Opt.scalar_struct);
tabappend = string(tabappend(:));
tabappend = strjoin(tabappend, "_");

tablefolder = figuredefine("tables");
t = table.analyses.eventuv(out, [], [], [], Opt.scalar_struct);
tmp = Option.patternNames(t.patterns);
t.patternNames = tmp(:);
assert(numel(unique(t.uv_components)) > 1, "Only one uv component found. Check your data.")
parquetwrite(fullfile(tablefolder, "eventuv" + tabappend + ".parquet"), t);
