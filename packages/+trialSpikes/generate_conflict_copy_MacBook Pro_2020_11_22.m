function [spikeSampleMatrix, spikeSampleTensor] = generate(spikeCountMatrix, timeAxis, Htrials, samplesPerTrial, numResult)
%
% Input
% -----
% spikes : numeric, time x neurons
%   spikeCountMatrix or spikeRateMatrix
%
% timeAxis: 1 x  time
%   Times correpsonding to the the time bins of spikes
%
% Htrials : numeric, trials  x 2
%   Start and stop time of each trial
%
% samplesPerTrial : scalar
%   Number of times to sample equidistant per trial
%
% Output
% ------
% spikeSampleMatrix : {neuron x (time*trial) per network pattern}
% spikeSampleTensor : {neuron x time x trial per network pattern}

% samplesPerTrial = 6;

[nNeurons,~] = size(spikeCountMatrix);

spikeSampleTensor = cell(1,numResult);
spikeSampleMatrix = cell(1,numResult);

for iPattern = progress(1:numResult, 'Title', 'Patterns')
    [nTrials,~] = size(Htrials{iPattern});
    currTensor = zeros(nNeurons,samplesPerTrial,nTrials);
    
    currPattern   = Htrials{iPattern};
    windows = length(currPattern(:,1));
    
    tq = zeros(1,samplesPerTrial);
    for i = progress(1:windows, 'Title', 'Times')
        currstart = currPattern(i,1);
        currstop  = currPattern(i,2);
        
        logicalIndexes = timeAxis>currstart & timeAxis<=currstop;
        
%         if sum(logicalIndexes) < 2
%             continue;
%         end
%                 if sum(logicalIndexes == 0)
%                     continue
%                 end
        
        t = timeAxis(logicalIndexes);
        %t_index = find(logicalIndexes);
        
        interval = (currstop-currstart)/samplesPerTrial;
        
        %Generte query points
        tq(1) = currstart+interval/2;
        for n = 2:samplesPerTrial
            tq(n) = tq(n-1)+interval;
        end
        %tq = linspace(currstart, currstop, samplesPerTrial+1);
        %tq = tq(1:end-1);
        
        % RY : added test assertions
        assert( min(tq) >= currstart, ...
            'violation: query times below window start')
        assert( max(tq) <= currstop, ...
            'violation: query times above window stop')
        
        % Per neuron interpolate
        for j = 1:nNeurons
            try
                neuron_spiking = spikeCountMatrix(j,logicalIndexes);
            catch
                keyboard
            end
            try
                
                interped_spiking = interp1(t,neuron_spiking, tq, 'linear');
                
            catch
                keyboard
                continue
            end
            currTensor(j,:,i) = interped_spiking;
        end
    end
    spikeSampleTensor{iPattern} = currTensor;
    temp = reshape(currTensor,[nNeurons,samplesPerTrial*nTrials]);
    temp(isnan(temp)) = 0;
    spikeSampleMatrix{iPattern} = temp;
end


end
