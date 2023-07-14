function out = match(Patterns, Option, r, varargin)
% MATCH Match source and target patterns with rank regressor
% or canonical correlation analysis (CCA) components.
%
%   out = match(Patterns, Option, r)
%
%   INPUTS
%       Patterns - struct with fields:
%           X_source     - source patterns
%           X_target     - target patterns
%           X_time       - time vector
%           index_source - index of source patterns
%           index_target - index of target patterns
%           rankRegress  - rank regressor analysis struct
%           cca          - canonical correlation analysis struct
%       Option - struct with fields:
%           method           - 'prod' or 'concat'
%           component_method - 'rank' or 'cca'
%           n_components     - number of components to use
%           smoothing        - samples or param to smooth()
%           verbose          - print out info
%    OUTPUTS
%       out - struct with fields:
%           activities       - matched activities
%           smooth_activities - smoothed activities



ip = inputParser();
ip.addParameter('method', 'prod', @(x) ismember(x, {'prod', 'concat'}));
ip.addParameter('component_method', 'rankRegress', @(x) ismember(x, {'rankRegress', 'cca'}));
ip.addParameter('n_components', 3, @isscalar);
ip.addParameter('smoothing', 400); % samples or param to smooth()
ip.addParameter('verbose', false, @islogical);
ip.parse(varargin{:});
Opt = ip.Results;

N = Opt.n_components;
out.method = Opt.method;
out.component_method = Opt.component_method;
out.activities = [];
out.smooth_activities = [];
out.time = [];

% source = Patterns.X_source;
% target = Patterns.X_target;
% time   = Patterns.X_time;
source_index = Patterns.index_source;
target_index = Patterns.index_target;
source = r.spikeCountMatrix(source_index,:);
target = r.spikeCountMatrix(target_index,:);
time = r.timeBinMidPoints;

% Pull out rank regressor
if strcmp(Opt.component_method, 'rankRegress')
    rr = Patterns.rankRegress;
    B_ = rr.B_;
    if isempty(B_)
        warning('No rank regressor found');
        return % no rank regressor
    end
    [u,~,v] = svd(B_);
    % sdiag = diag(s);
elseif strcmp(Opt.component_method, 'cca')
    cca = Patterns.cca;
    if isempty(cca)
        warning('No cca found');
        return % no cca
    end
    u = cca.u;
    v = cca.v;
else
    error('Unknown method');
end
if Opt.verbose
    disp("Size source: " + size(source))
    disp("Size u: " + size(u))
    disp("Size target: " + size(target))
    disp("Size v: " + size(v))
end

% Project onto rank regressor
if strcmp(Opt.method, 'concat')
    both = [source; target];
    % TODO: center both?
    match_inputoutput = [u(:,1:N); v(:,1:N)];
    activities = match_inputoutput' * (both);
    activities = activities - mean(activities,2);
    activities = abs(activities);
    out.activities = activities;
elseif strcmp(Opt.method, 'prod')
    % Closer to what occurs in SVD and regression
    % requires product of activations > 0 for match
    activation_source = u(:,1:N)'*source;
    activation_target = v(:,1:N)'*target;
    activities = activation_source .* activation_target;
    out.activities = activities;
    out.activaties_source = activation_source;
    out.activaties_target = activation_target;
else
    error('Unknown method');
end

% Smooth
smooth_activities = zeros(size(activities));
for i = progress(1:size(activities,1))
    smooth_activities(i,:) = smooth(activities(i,:), Opt.smoothing);
end
out.smooth_activities = smooth_activities;
out.time = time;
