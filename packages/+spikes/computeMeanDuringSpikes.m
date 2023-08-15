function means = computeMeanDuringSpikes(times_spiking, behavior, column_name, varargin)
% COMPUTEMEDIANDURINGSPIKES Compute median value of a behavior column during spikes
%

    ip = inputParser();
    ip.addParameter('quantile_vec', [], @isvector);
    ip.addParameter('quantile', 0.01, @isscalar);
    ip.parse(varargin{:});
    Opt = ip.Results;

    % Check if column_name exists in behavior table
    if ~ismember(column_name, behavior.Properties.VariableNames)
        error('Column name does not exist in behavior table.');
    end

    numNeurons = numel(times_spiking);
    means = zeros(numNeurons, 1);

    for n = progress(1:numNeurons, 'Title', 'Computing median values')

        % Extract spike times for the current neuron
        neuronSpikeTimes = times_spiking{n};

        % Extract values of the behavior column at the times of spikes
        % [~, ia] = ismembertol(behavior.times, neuronSpikeTimes, 1e-10, 'DataScale', 1);
        behavior_inds = interp1(behavior.time, 1:numel(behavior.time), neuronSpikeTimes, 'nearest', 'extrap');
        inds = abs(behavior.time(behavior_inds) - neuronSpikeTimes(:)) > 0.033;
        behavior_inds(inds) = [];
        spikeBehaviorValues = behavior.(column_name)(behavior_inds);
        if isempty(spikeBehaviorValues)
            medians(n) = NaN;
            continue
        end

        % Quantile filter
        if isempty(Opt.quantile_vec)
            Opt.quantile_vec = spikeBehaviorValues(:);
        end
        Q1 = quantile(Opt.quantile_vec, Opt.quantile);
        Q2 = quantile(Opt.quantile_vec, 1-Opt.quantile);
        spikeBehaviorValues(Opt.quantile_vec < Q1) = [];
        spikeBehaviorValues(Opt.quantile_vec > Q2) = [];
        % Compute the median value
        means(n) = mean(spikeBehaviorValues);

    end
end
