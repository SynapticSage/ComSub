function Out = timeVarying(Patterns, Option, r, varargin)
% Components = timeVarying(Patterns, Option)
%
% This function takes the patterns and performs the time varying analysis
% on them. It returns a struct with the critical components for each
% behavior and the subspaces for each behavior.
%
% TODO: 
% 1. How much each cell matches the subspace
% 2. How much all match
%    - old RRR method
%    - new RRR method
%    - CCA method
% 3. Relationship to various behavioral variables
%   - speed
%   - acceleration
%   - choice point (out in)
%   - reward location (out in)
%
% Inputs:
% -----
%   Patterns: struct with the patterns for each partition
%   Option: struct with the options for the analysis
%
% Outputs:
%   Components: struct with the critical components for each behavior and
%   the subspaces for each behavior

ip = inputParser();
ip.addParameter('method', 'prod');
ip.addParameter('componentMethod', 'rankRegress');
ip.parse(varargin{:});
Opt = ip.Results;

disp("Running time varying analysis")
tic
Const = option.constants();

if ~isfield(Patterns, 'rankRegress') || isempty(Patterns(1).rankRegress)
    error("Patterns must have the rankRegress field")
end

if Opt.componentMethod == "rankRegress"
    outkey = "rrr";
elseif Opt.componentMethod == "cca"
    outkey = "cca";
else
    error("Unknown component method")
end


clear Components
% Components = repmat(Components, ...
%                     [Option.numPartition, Option.waysOfPartitions]);

for i = progress(1:numel(Patterns), 'Title', 'Time varying analysis')

    P = Patterns(i); 
    if i > 1
        O = Out(i);
    else
        O = struct();
    end
    O.directionality = P.directionality;
    tmp = split(P.directionality, '-');
    O.source = tmp{1};
    O.target = tmp{2};

    % old method
    % ----------------
    % old = analysis.temporal.oldMatching(p, Option, r, target, animal_behavior,...
    %                               unique_times, throwout_times);

    % new method
    % ----------------
    O.(outkey) = analysis.temporal.match(P, ...
                Option, r,...
                'method', Opt.method,...
                'component_method', Opt.componentMethod...
                );
    O.(outkey).method = Opt.method;
        % O.("cca_"+method{1}) = analysis.temporal.match(P, ...
        %             Option, r,...
        %             'method', method{1},...
        %             'component_method', 'cca'...
        %             );
    O.name = P.name;

    if i == 1
        disp("Initializing output struct")
        Out = repmat(nd.emptyLike(O), size(Patterns));
    end
    Out(i) = O;
end

disp("Finished time varying analysis in " + string(toc) + " seconds")
