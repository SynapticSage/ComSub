function [results] = core(x_area1, x_area2, varargin)
% core canonical correlation analysis
%   [results] = core(x_area1, x_area2) performs canonical correlation
%   analysis on the two matrices x_area1 and x_area2. the matrices should
%   have the same number of rows (observations) and different number of
%   columns (variables). the function returns a structure with the
%   following fields:
%       a: canonical coefficients for area1
%       b: canonical coefficients for area2
%       r: canonical correlations
%       u: canonical scores for area1
%       v: canonical scores for area2
%       stats: statistics for the test of h0: r = 0
%       area1_scores: top canonical scores for area1
%       area2_scores: top canonical scores for area2
%       cross_area_corr: correlation between the top canonical scores
%

% ip = inputParser;
% ip.addParameter('ploton', false, @islogical);
% ip.parse(varargin{:});
% opt = ip.Results;
d=@double;
usecv = true;

% run cca
if ~usecv
    [a, b, r, u, v, stats] = canoncorr(d(x_area1)', d(x_area2)');
    R = nan;
else
    out = analysis.cca.cvcca(d(x_area1)', d(x_area2)');
    u = out.U;
    v = out.V;
    a = out.A;
    b = out.B;
    stats = out.pval;
    R = out.val;
    r = diag(R);
end

% correlation between the first scores
cross_area_corr = corr(u(:,1), v(:,1));

% top cv scores
area1_scores = u(:,1);
area2_scores = v(:,1);

results = struct(...
    'a', a, ...
    'b', b, ...
    'r', r, ...
    'R', R, ...
    'u', u, ...
    'v', v, ...
    'stats', stats, ...
    'area1_scores', area1_scores, ...
    'area2_scores', area2_scores, ...
    'cross_area_corr', cross_area_corr ...
);

% make the plotf
scatter(area1_scores, area2_scores); % scatter plot of the scores
hold on; % keep the scatter plot when adding the next plot
plot([-1, 1], [-1, 1], 'k--'); % add the x=y line (adjust the range as needed)
xlabel('area1 cross-activity');
ylabel('area2 cross-activity');
title('cross-activity');

end
