function medians = computeMedianDuringSpikes(times_spiking, behavior, column_name)
    % Check if column_name exists in behavior table
    if ~ismember(column_name, behavior.Properties.VariableNames)
        error('Column name does not exist in behavior table.');
    end

    numNeurons = numel(times_spiking);
    medians = zeros(numNeurons, 1);

    for n = 1:numNeurons
        % Extract spike times for the current neuron
        neuronSpikeTimes = times_spiking{n};

        % Extract values of the behavior column at the times of spikes
        % [~, ia] = ismembertol(behavior.times, neuronSpikeTimes, 1e-10, 'DataScale', 1);
        behavior_inds = interp1(behavior.time, 1:numel(behavior.time), neuronSpikeTimes, 'nearest', 'extrap');
        inds = abs(behavior.time(behavior_inds) - neuronSpikeTimes) > 0.033;
        behavior_inds(inds) = [];
        spikeBehaviorValues = behavior(behavior_inds, column_name);

        % Compute the median value
        medians(n) = median(spikeBehaviorValues, 'omitnan');
    end
end
