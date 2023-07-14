function [result]  = optimizeOptionsScript(Patterns, Raw, Options)
% this script evaluates how ideal an option struct is 
% by variable result = corr (theta - speed) + corr(theta - ripple)
%                      - avgperf(hpc-pfc)
% avgperf is the averge performance of hpc predicting pfc across the three
% pattern rhythms



corrMatrix = params.getPatternSpeedCorr(Option.animal,...
                                        Raw.Events.H,...
                                        Raw.Events.Hvals,...
                                        Raw.Events.Htimes);

theta_sp     = corrMatrix(2,1);
theta_ripple = corrMatrix(4,2);


% Distribution of performance across partitions
[full_model_performance, dist] = params.getPredictionPerf(Patterns, nTarget, nSource, Option.numPartition);

% Visualize
fig('params.getPredictionPerf')
histogram(dst)
perf = sum(full_model_performance(2:2:end))/3;

% Compute objective
optimizationResult = theta_sp + theta_ripple - perf;


result.theta_sp = theta_sp;
result.theta_ripple = theta_ripple;
result.full_model_performance = full_model_performance;
result.dist = dist;
result.optimizationResult = optimizationResult;
