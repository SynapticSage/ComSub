function out = getSpikeTrain(animal, timebinSize, samplingRate)
% GETSPIKES Summary of this function goes here
%
% Input
% -----
% animal: name of the animal
%
% timebinSize : float
%   Size of time bins in seconds for spikeCountMatrix and
%   spikeRateMatrix. if given this option alone this determines both
%   sampling rate and window size.
%
% samplingRate : float
%   if given this (in addition to timebinSize), a time bin size window
%   is swept along the time axis collecting samples at the
%   samplingRate.
%
%  if not given, the timebinSize is used as the sampling rate
%   and the window size is set to 1 second.
%
% Output
% -----
% timeBinStartEnd: axis of start and end times of each time bin
%
% timeBinMidPoints: an array of all the time bin's midpoints
%
% times_spiking : cell
%   For every neuron, one list of spike times across the sessions,
%   stored in a cell
%
%
%
% spikeCountMatrix : timebins x neurons matrix
%   Spike counts per timebin
%
% spikeRateMatrix : timebins x neurons matrix
%   Spike rate per timebin
%
% areaPerNeuron : 1 x neurons string "" matrix
%   List of area per neuron. Used to later separate the columns
%   corresponding to each brain area. This function gathers that information
%   from the cellinfo file associated with the filename, by exacting animal
%   name from filename and appending 'cellinfo.mat'
%
% cell_index: the tetrode*cell tuple of all cells in the interactions
%end
%
%% setting up inputs
%     filename = "JS15spikes01.mat";
%     timebinSize = 0.1;

disp("...loading spikes for " + animal)
tic

load(animal + "spikes01.mat");
% % In this case using for-looping integers is easier than raw forlooping the
% % underlying structures. Easier than forlooping the underlying structures


%------------------------------------------------------------------
% to be implemented:
% given an animal get to the directory and get the spikes structure
%     animal = append(extractBetween(filename, 1, 4),"cellinfo.mat");
load (animal + "cellinfo.mat");

%------------------------------------------------------------------

%% select unique neurons
session = 1;
spikes = spikes{session}; % single day epoches
cellinfo = cellinfo{session};
num_tetrode = numel(spikes{1});

cell_index = [];

for epoch = 1:numel(spikes)
    epochSpikes = spikes{epoch};
    for tetrode = 1:numel(epochSpikes)
        tetrodeEpSpikes = epochSpikes{tetrode};
        for neuron = 1:numel(tetrodeEpSpikes)
            neuronTetEpoch = tetrodeEpSpikes{neuron};
            if ~isempty(neuronTetEpoch) && isfield(neuronTetEpoch,'data')...
                    && ~isempty(neuronTetEpoch.data)
                cell_index(end+1,:) = [tetrode,neuron];
            end
        end
    end
end

% select unique cells
cell_index = unique(cell_index(:,[1 2]),"rows");

size_cellindex = size(cell_index);
num_cells = size_cellindex(1);


areaPerNeuron = string.empty;
% time_spiking should contain the unique cells, why the loop does not work? 
times_spiking = cell(1,num_cells);

%% get information about each neuron
for epoch = 1:numel(spikes)
    for tetrode = 1:numel(spikes{epoch})
        tetrodeEpSpikes = spikes{epoch}{tetrode};
        for neuron = 1:numel(tetrodeEpSpikes)
            if (numel(spikes{epoch}{tetrode}{neuron})~=0)
                for i = 1:num_cells
                  %  try
                        if neuron == cell_index(i,2) && tetrode == cell_index (i,1)
                            if isempty(spikes{epoch}{tetrode}{neuron}.data)
                                continue
                            else
                                % disp("...getting spike times for neuron " + neuron + " on tetrode " + tetrode + " on epoch " + epoch)
                                neuron_data = spikes{epoch}{tetrode}{neuron}.data(:,1);
                                times_spiking{i} = [times_spiking{i}, neuron_data'];
                                
                                areaPerNeuron(i) = cellinfo{epoch}{tetrode}{neuron}.area;
                            end
                        end
%                     catch
%                         keyboard
%                     end
                end
            end
        end
    end
end

%% generate spike count/rate matrices
start_time = intmax;
end_time = -1;

% find the start & end time of the day
disp("...finding start and end times of the day")
for i = 1:num_cells
    if isempty(times_spiking{i})
        continue
    end
    
    try
        curr_start = times_spiking{i}(1);
        curr_end = times_spiking{i}(end);
    catch
        keyboard
    end
    
    if curr_start<start_time
        start_time = curr_start;
    end
    
    if curr_end>end_time
        end_time = curr_end;
    end
end

%%
if nargin == 3 && ~isempty(samplingRate) %RY added new option for sampling rate + window size, copy-pasting some of your code and modifying
    
    samplingPeriod = 1/samplingRate;
    timeBinStartEnd = start_time:samplingPeriod:end_time;
    nTimes = length(timeBinStartEnd);
    spikeCountMatrix = zeros(num_cells,nTimes);
    
    % chunk size is the number of time bins to process at a time -- purely for
    % the amount of time to process each loop
    chunkSize = 500;
    
    %         spikeCountSlice = zeros(1, nTime);
    %         for i = progress(1:num_cells,'Title','cells') % surrounding a 1:num_cells with progress() adds a progress bar
    %             for t = progress(1:chunkSize:numel(timeBinStartEnd),'Title','TimeChunks')
    %                 windowStarts = (timeBinStartEnd(t:min(t+chunkSize, nTime)) - timebinSize/2);
    %                 windowStops =  (timeBinStartEnd(t:min(t+chunkSize, nTime)) + timebinSize/2);
    %                 spikeCountSlice(t:min(t+chunkSize,nTime)) = ...
    %                     sum( times_spiking{i}' >= windowStarts & ...
    %                     times_spiking{i}' < windowStops);
    %             end
    %             spikeCountMatrix(i, :) = spikeCountSlice;
    %         end
    %         spikeRateMatrix = spikeCountMatrix/timebinSize;

    % windowStarts/Stops determines the start and stop times of each time bin
    windowStarts = (timeBinStartEnd - timebinSize/2);
    windowStops =  (timeBinStartEnd + timebinSize/2);
    p = ProgressBar(num_cells, ...
        'Title', 'cells');
    %'IsParallel',true);
    p.setup([],[],[]);
    cleanupFunction = onCleanup(@(x) ProgressBar.deleteAllTimers());
    for i = 1:num_cells % surrounding a 1:num_cells with progress() adds a progress bar
        one_cell_spiking = times_spiking{i};
        nSpikes = numel(one_cell_spiking);
        spikeCountSlice = zeros(1, nTimes);
        for t = progress(1:chunkSize:numel(times_spiking{i}) ,'Title','SpikeChunks')
            spikeChunk = one_cell_spiking(t:min(t+chunkSize, nSpikes));
            spikeCountSlice = spikeCountSlice + ...
                sum( spikeChunk' >= windowStarts & ...
                spikeChunk' < windowStops);
        end
        spikeCountMatrix(i, :) = spikeCountSlice;
        p.step([], [], []);
    end
    spikeRateMatrix = spikeCountMatrix/timebinSize;
    
else % standard option, just window size alone
    
    timeBinStartEnd = start_time:timebinSize:end_time;
    timeBinMidPoints = zeros(1,numel(timeBinStartEnd)-1);
    for i = 1:numel(timeBinStartEnd)-1
        timeBinMidPoints(i) = (timeBinStartEnd(i) + timeBinStartEnd(i+1))/2;
    end
    
    spikeCountMatrix = zeros(num_cells,length(timeBinStartEnd)-1);
    spikeRateMatrix = zeros(num_cells,length(timeBinStartEnd)-1);
    
    for i = 1:num_cells
        [spike_count,timeBinStartEnd] = histcounts(times_spiking{i},timeBinStartEnd);
        spikeCountMatrix(i,:) = spike_count;
        spikeRateMatrix(i,:) = spike_count/timebinSize;
    end
end

%% get the session type per bin
% whether each bin belongs to sleep (0) or run (1)

startStop = getEpochTimes(animal);
sessionTypePerBin = zeros(1,numel(timeBinMidPoints));
epochPerBin = zeros(1,numel(timeBinMidPoints));
%%
for i = 1:numel(timeBinMidPoints)
    currTime = timeBinMidPoints(i);
    for j = 1:size(startStop,2)
        if currTime >= startStop(1,j) && currTime < startStop(2,j)
            if mod(j,2) == 0 
                sessionTypePerBin(i) = 1;
            else
                sessionTypePerBin(i) = 0;
            end
            epochPerBin(i) = j;
        end
    end
end

out.spikeCountMatrix = spikeCountMatrix;
out.spikeRateMatrix = spikeRateMatrix;
out.timeBinStartEnd = timeBinStartEnd;
out.timeBinMidPoints = timeBinMidPoints;
out.times_spiking = times_spiking;
out.areaPerNeuron = areaPerNeuron;
out.cell_index = cell_index;
out.sessionTypePerBin = sessionTypePerBin;
out.epochPerBin = epochPerBin;

disp("...done " + toc + " seconds")
