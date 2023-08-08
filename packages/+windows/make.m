function [cellOfWindows, cutoff, opposite_cutoff] = ...
    make(times, threshold, H, Option, varargin)
% H: m*n matrix - m is number times taken (aligned to times, same
%    length as times) and n is the number of rhythms we are interested in
%
% times: time points where the powers of different rhythms were taken
%
% threshold: percentile for which the powers are chosen
%
% Option: Option structure, see option.defualts() for details
%
% if Opt.quantile given and H, then the Opt.quantile sequencec is used to find
% the quantile() and H is the sequence that the cutoff is applied to in order
% to attain windows.
%
% 6/23 - RY - added the ability to throw out outlier values
% 6/23 - RY - vectorized the binary network pattern detection (ensured it
%             outputs the same pattern as before)
% 6/23 - RY - added the ability to check for positive derivative within the
%             window -- useful for theta and delta, ensuring majority of the
%             window does not contain pattern diving below the threshold
%             (or if higherThanQuantile is false, above the threshold)

winSize        = Option.winSize;
posDeriviative = Option.positiveDerivativeCheck;

% Optional argumuents
% -------------------
ip = inputParser;
ip.addParameter('intraWindow_binCenters', []);       % How much to shift between samples .. if not given, it's the winSize
ip.addParameter('intraWindow_binCenter_edges', []);         % How much to shift between samples .. if not given, it's the winSize
ip.addParameter('thresholdMethod','quantile'); % how do we interpret the numbers given in threshold, 'quantile' as a quantilee, and 'raw' as a raw threshold
ip.addParameter('quantile',[]);                % RY if not empty, then it calculates quantile of this variable instead of for H
ip.addParameter('higherThanQuantile', true);   % when there is a single threshold value, this encodes the directionality of the test, >  if true, < if false
ip.addParameter('outlierQuantile', [0, 1]);% quantile to use for outlier detection
ip.addParameter('positiveDerivativeCheck', false) % if not empty, then it checks if the derivative is positive for some range within the window
ip.parse(varargin{:})
Opt = ip.Results;

posDerivativeMode = ~isempty(posDeriviative) && ~isnan(any(posDeriviative)) ...
    && Opt.positiveDerivativeCheck;

boolcheck = [isempty(Opt.intraWindow_binCenters), ...
             isempty(Opt.intraWindow_binCenter_edges)];
assert(all(boolcheck==0) || all(boolcheck==1),...
    'Either provide both intraWindow_binCenters and intraWindow_binCenter_edges or none');

%% 1. calculate the quantiles
%% --------------------------
[~,nPatterns] = size(H);
cutoff = zeros(numel(threshold),nPatterns);
opposite_cutoff = zeros(numel(threshold),nPatterns);

% Compute quantile cutoffs either from H or another given variable
for i = 1:nPatterns
     switch Opt.thresholdMethod
        case 'quantile'
            if ~isempty(Opt.quantile) % USER DEFINED QUANTILE VECTOR
                cutoff(:,i) = quantile(Opt.quantile(:,i),threshold);
                opposite_cutoff(:,i) = quantile(Opt.quantile(:,i),1-threshold);
            else % QUANTILE OF H
                cutoff(:,i) = quantile(H(:,i),threshold);
                opposite_cutoff(:,i) = quantile(H(:,i),1-threshold);
            end
        case 'raw'
            cutoff(:,i) = threshold;
            opposite_cutoff(:,i) = nan;
        otherwise
            error('Invalid method')
    end
end

upperOutlierCutoff = quantile(H, Opt.outlierQuantile(2));
lowerOutlierCutoff = quantile(H, Opt.outlierQuantile(1));
upperOutlierCutoffRepeat = repmat(upperOutlierCutoff, [size(H, 1), 1]);
lowerOutlierCutoffRepeat = repmat(lowerOutlierCutoff, [size(H, 1), 1]);

%% --------------------------------
%% 2. Apply quantiles to the powers
%% --------------------------------
%marking the H_matrix with 1 if the power is above quantile and 0 if
LOW = 1;
HIGH = 2;
is_cutoff_range = size(cutoff, 1) == 2;
assert(size(cutoff,1) == 1, "Cutoff must be a range (2 numbers) or a min or max (1 number)")
% H_isPattern = false(size(H));

% Repeat the cutoff array to match size of H
cutoffRepeat = repmat(cutoff, [size(H, 1), 1]);
cutoffOppositeRepeat = repmat(opposite_cutoff, [size(H, 1), 1]);

if is_cutoff_range
    disp('selecting range')
    H_isPattern = (H >= cutoffRepeat(:,LOW)) & (H < cutoffRepeat(:,HIGH));
else
    if Opt.higherThanQuantile == 0 || Opt.higherThanQuantile == false
        disp("selecting low " + threshold)
        H_isPattern = H < cutoffRepeat;
    elseif Opt.higherThanQuantile == 1 || Opt.higherThanQuantile == true
        disp("selecting high = " + threshold)
        H_isPattern = H >= cutoffRepeat;
    elseif Opt.higherThanQuantile == -1
        disp("selecting middle " + threshold + " " + opposite_threshold)
        if opposite_cutoff(1) < cutoff(1)
            H_isPattern = (H >= cutoffOppositeRepeat) & (H < cutoffRepeat);
        else
            H_isPattern = (H >= cutoffRepeat) & (H < cutoffOppositeRepeat);
        end
    else
        error('Invalid higherThanQuantile value')
    end
end

% Marking as false for H values outside upper and lower quantile range
H_isPattern = H_isPattern & (H > lowerOutlierCutoffRepeat) &...
                            (H < upperOutlierCutoffRepeat);

if isscalar(winSize)
    winSize = [-winSize, winSize];
elseif iscell(winSize)
    winSize = [winSize{:}];
end


%% -------------------------
%% 3. map to the time points
%% -------------------------
cellOfWindows = cell(1, length(cutoff));
% the time points we want, 1st column starting and 2nd ending
for pattern = 1:length(cutoff)

    switches = diff(H_isPattern(:,pattern));
    switches = (switches == 1);
    temp = find(switches);
    temp = (switches == 1);
    if iscolumn(times)
        times = times';
    end
    
    trigger_times = times(temp)'; % Capture triggers
    starts = trigger_times+winSize(1);
    stops  = trigger_times+winSize(2);

    % if we want to check for positive derivative (i.e. the derivative is
    % positive for some range within the window)
    if posDerivativeMode 
        % determine start and stop times to check for positive derivative
        % (i.e. the derivative is positive between start and stop)
        derstarts = trigger_times+posDeriviative(1);
        derstops  = trigger_times+posDeriviative(2);
        % find the indices of the time indices for ranges to check
        derstarts = interp1(times, 1:length(times), derstarts, 'nearest');
        derstops  = interp1(times, 1:length(times), derstops,  'nearest');
        derranges = [derstarts, derstops];
        derranges = derranges(all(~isnan(derranges),2),:);
        correct_direction = false(size(derranges,1),1);
        for i = 1:size(derranges,1)
            % check if cumulative direction is positive 
            if is_cutoff_range || Opt.higherThanQuantile == true || Opt.higherThanQuantile == 1
                if H(derranges(i,2),pattern) - H(derranges(i,1),pattern) > 0
                    correct_direction(i) = true;
                end
            elseif Opt.higherThanQuantile == false || Opt.higherThanQuantile == 0
                if H(derranges(i,2),pattern) - H(derranges(i,1),pattern) < 0
                    correct_direction(i) = true;
                end
            else % more middle range windows, we're neither looking for positive or negative derivative
                correct_direction = true(size(derranges,1),1);
            end
        end
        disp("Derivative mode for Pattern " + pattern)
        disp("...fraction of positive derivative windows: " + mean(correct_direction))
    else
        correct_direction = true(size(starts));
    end

    if isempty(Opt.intraWindow_binCenters) % Cell of windows is windows x 2 (starts and ends of whole window)
        windows = [starts, stops];
    else % Cell of windows is windows x 2 x timeChunks (starts and ends of chunks)
        shifts = 0 : Opt.intraWindow_binCenters : sum(abs(winSize));
        chunkCenters = starts + shifts;
        windows = chunkCenters + Opt.intraWindow_binCenter_edges;
    end
    cellOfWindows{pattern} = windows(correct_direction,:);

end
