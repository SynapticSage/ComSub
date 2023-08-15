function [bestloc] = computeOccupancyNormalizedMean(times_spiking, behavior, column_name, varargin)
    ip = inputParser();
    ip.addParameter('quantile_vec', [], @isvector);
    ip.addParameter('quantile', 0.01, @isscalar);
    ip.parse(varargin{:});
    % Opt = ip.Results;

    % Check if column_name exists in behavior table
    if ~ismember(column_name, behavior.Properties.VariableNames)
        error('Column name does not exist in behavior table.');
    end

    numNeurons = numel(times_spiking);
    bestloc = zeros(numNeurons, 1);

    % Calculate total time and time for each behavior bin
    dtime = diff(behavior.time);
    uniqueBins = unique(behavior.(column_name));
    totalTime  = sum(dtime(dtime > 0 & dtime < 1));
    timePerBin = arrayfun(@(x) sum(behavior.time(behavior.(column_name) == x)), uniqueBins);
    
    for n = progress(1:numNeurons, 'Title', 'Computing occupancy normalized mean values')
        % Extract spike times for the current neuron
        neuronSpikeTimes = times_spiking{n};
        
        % Extract bins of the behavior column at the times of spikes
        behavior_inds = interp1(behavior.time, 1:numel(behavior.time), neuronSpikeTimes, 'nearest', 'extrap');
        spikeBehaviorBins = behavior.(column_name)(behavior_inds);
        
        spikeCountsPerBin = histc(spikeBehaviorBins, uniqueBins);
        firingRates = spikeCountsPerBin ./ timePerBin;  % Firing rate per bin
        firingRates = firingRates * 100;
        
        % Normalize firing rates using the occupancy normalization formula
        normalizedFiringRates = firingRates .* (totalTime ./ timePerBin);
        normalizedFiringRates = smooth(normalizedFiringRates, 3);

        [~,bestloc(n)] = max(normalizedFiringRates);
    end
end

