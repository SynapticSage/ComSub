function [Patterns] = ccaAnalysis(Patterns, Option)
% ccaAnalysis - Apply CCA to each pattern struct
%   Patterns = ccaAnalysis(Patterns)
%
% INPUTS:
%   Patterns - [1 x 1 struct] see partitionPatterns.m for fields
%
% OUTPUTS:
%   Patterns - [1 x 1 struct] with additional fields:
%       cca - [1 x 1 struct] see cca.core for fields

disp('Running CCA analysis...');
tic;

intra_only_overall = true;

% Iterate over each pattern and apply CCA
for n = 1:numel(Patterns)
    p = Patterns(n);

    curr_area1 = p.X_source;
    curr_area2 = p.X_target;
    equals = isequal(p.X_source, p.X_target);

    disp("CCA: " + p.name + " isequal=" + equals);
    if intra_only_overall  &&equals  &&p.name  ~= "Overall"
        disp("Skipping...")
        continue;
    else
        disp("Running...")
    end


    if(mean(curr_area1, 'all') < 100*eps())
        warning('centering...');
        curr_area1 = curr_area1 - mean(curr_area1, 2);
        curr_area2 = curr_area2 - mean(curr_area2, 2);
    end
    
    % Assuming cca.core function takes two arguments: X_area1 and X_area2
    result = analysis.cca.core(curr_area1, curr_area2, 'area1_eq_area2', equals);
    if equals  && p.name  =="Overall"
        disp("target 2nd dim size " + size(Patterns(end).X_target, 1))
        result2 = analysis.cca.core(p(end).X_target, Patterns(end).X_target, 'area1_eq_area2', equals);
        result.b = result2.b;
        result.v = result2.v;
    end
    p.cca = result;
    
    % Assign the updated struct back to the original array
    Patterns = nd.set(Patterns, n, p);
end

disp(['CCA analysis completed in ' num2str(toc) ' s']);
