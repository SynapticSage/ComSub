function Events = ThetaDeltaRipple(Events, Option)
% [cellOfWindows, cutoffs] = ThetaDeltaRipple(Events, Option)
%   Generates windows for theta, delta, and ripple events
%   INPUTS:
%       Events: struct with fields:
%           times: 1 x nEvents vector of event times
%           H: nEvents x nPatterns matrix of H values
%           Hvals: nEvents x nPatterns matrix of H values
%       Option: struct with fields:
%           patternNames: 1 x nPatterns cell array of pattern names
%               (e.g. {'SWR', 'Ripple', 'Theta', 'Delta'})
%           winSize: size of window in seconds
%               (e.g. 0.5)
%           quantileToMakeWindows: quantile of H values to make windows
%               (e.g. 0.9)
%           equalWindowsAcrossPatterns: true/false
%               meaning if true, equalizes number of windows across patterns
%           generateH: cell array of strings of H to generate
%               meaning if {'fromCoherence', 'fromWpli'}, then will generate
%               H from coherence and wpli
%           lowerControl: true/false
%               meaning if true, will generate control windows from lower
%               quantile of H values
%           oldControlBehavior: true/false
%               meaning if true, will generate control windows from lower
%               quantile of H values
%           singleControl: true/false
%               meaning if true, will generate control windows from lower
%               quantile of H values
%   OUTPUTS:
%       cellOfWindows: 1 x nPatterns cell array of windows
%       cutoffs: nPatterns x 1 vector of cutoffs

disp("Generating windows for events")
tic;

const = option.constants();

THETA  = const.THETA;
DELTA  = const.DELTA;
RIPPLE = const.RIPPLE;

nPatterns = numel(Option.patternNames);

%%--------------------
%% 1s: THETA AND DELTA
%%--------------------
[cellOfWindows, cutoffs] = windows.make(Events.times, ...
    Option.quantileToMakeWindows, Events.H(:,THETA:DELTA), Option,...
    'positiveDerivativeCheck', true, ...
    'outlierQuantile', Option.thetadelta_outlierQuantile);

%%----------------
%% 1b :RIPPLES
%%----------------
coherent_patterns = any(contains(Option.generateH, ["fromCoherence","fromWpli"]));
if coherent_patterns
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(...
        Events.times,...
        Option.quantileToMakeWindows, ...
        Events.Hvals(:,RIPPLE),... 2023, rather than rip times with high coh (H), just high ripple coherence (Hvals)
        Option,...
        'quantile', Events.Hvals(:,RIPPLE),'higherThanQuantile', true); 
    % RY: quantile needs to be hvals for ripple coherence/wpli threshold to
    % be correct, but timesd computed from Events.H such that non-ripple
    % times thrown out
    % Sample up to match lowest of theta/delta (ripple coherence will not be
    % compared to theta/delta; so let's not allow it to limit theta/delta
    % windows)
else
    % LITERALLY JUST RIP TIMES, NO THRESHOLD
    threshold = 1; % ripple is happening
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(...
        Events.times, ...
        threshold,...
        Events.H(:,RIPPLE),...
        Option, 'threshold', 'raw','higherThanQuantile', true);
end

windows.countMessage(cellOfWindows, Option.patternNames,...
    'message', 'initial window creation')

%%----------------------
%% 2: EQUALIZE N OF WINDOWS
%%----------------------
if Option.equalWindowsAcrossPatterns == true
    cellOfWindows = windows.equalizeWindowsAcrossPatterns(cellOfWindows);
end
windows.countMessage(cellOfWindows, Option.patternNames, ...
    'message', 'equalize windows')


numWindowsCut = size(cellOfWindows{1},1);
%windows.printWindowOverlap(cellOfWindows, Option.patternNames);

% -----------------------
% Controls: THETA AND DELTA
% -----------------------
if Option.lowerControl
    quantileControl = 1 - Option.quantileToMakeWindows;
else
    quantileControl = Option.quantileToMakeWindows;
end
if Option.oldControlBehavior
    Hc =  control.generatePatternShuffle(Events.H(:,1:3), Events.times, cellOfWindows); % add control patterns;
else
    Hc = Events.H;
end

[Hc_cellOfWindows, Hc_cutoffs] = ...
    windows.make(Events.times,  quantileControl,...  % add windows of control patterns
    Hc(:,THETA:DELTA), Option,... % Selects less than quantile
    'outlierQuantile', Option.thetadelta_outlierQuantile,...
    'positiveDerivativeCheck', true, ...
    ... under old controls, flag will select from shuffled times, under new, flag will select
    ... lower than quantile
    'higherThanQuantile', Option.oldControlBehavior);


% -----------------------
% Controls: RIPPLES
% -----------------------
coherent_patterns = any(contains(Option.generateH, ["fromCoherence","fromWpli"]));
if coherent_patterns
    [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(...
        Events.times,...
        quantileControl, ...
        Hc(:,RIPPLE), ...
        Option,...
        'quantile', Events.Hvals(:,RIPPLE),'higherThanQuantile', Option.oldControlBehavior);
else
    if Option.oldControlBehavior
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Events.times, ...
            1,   Hc(:,RIPPLE), Option, 'threshold', 'raw','higherThanQuantile', true); %RY quantile won't work because these are raw
    else
        % QUANTILE THRESH, PURE RIPPLE BAND, UNFILT BY RIP TIMES
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(...
            Events.times, ...
            quantileControl,...
            Events.Hvals(:,RIPPLE),...
            Option, ...
            'threshold', 'quantile','higherThanQuantile', Option.oldControlBehavior); 
    end
end

% ----------------------------------
% Clean, Remove overalapping windows
% ----------------------------------
% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:length(cellOfWindows)
    curr = windows.removeOverlapsBetweenPattern(...
        cell2mat(cellOfWindows(:,pattern)), cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end

% -----------------------------------
% Merge many controls into 1 control?
% -----------------------------------
% % Merge into one
control_start = length(cellOfWindows)+1;
control_end = length(cellOfWindows)+length(cellOfWindows);
cellOfWindows(control_start:control_end) = Hc_cellOfWindows;
cutoffs = [cutoffs,Hc_cutoffs];

% -------------------------------
% If any cellOfWindows is empty, throw an error
% -------------------------------
if any(cellfun(@isempty, cellOfWindows))
    error("A network pattern has no windows. Check your settings...")
end

% -----------------------------------------
% Ensure each pattern has equal # of window
% -----------------------------------------
% Equalize trials/windows for each pair of patttern-controlPattern
[cellOfWindows, warnedEmptyControls] =...
    control.equalizePatternControl(cellOfWindows);

% pick control pattern that actually contains controls, would break if all
% three are empty...
if Option.singleControl  &&warnedEmptyControls
    for iPossibleControl = nPatterns+1:nPatterns*2
        if ~isempty(cellOfWindows{iPossibleControl})
            cellOfWindows{nPatterns+1} = cellOfWindows{iPossibleControl};
            disp("here")
        end
    end
end

% -------------------------------
% If any cellOfWindows is empty, throw an error
% -------------------------------
if any(cellfun(@isempty, cellOfWindows))
    error("A network pattern has no windows. Check your settings...")
end

Events.cellOfWindows = cellOfWindows;
Events.cutoffs = cutoffs;
Events.nWindows      = cellfun(@(x) size(x, 1), cellOfWindows);

disp("")
disp("Windows generated " + num2str(numWindowsCut) + " windows cut")
disp("  time elapsed: " + num2str(toc) + " seconds")
disp("  cutoffs: " + num2str(cutoffs))
