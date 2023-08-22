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

intra_only_overall      = true;
remove_inter_from_intra = true;

% Iterate over each pattern and apply CCA
for n = 1:numel(Patterns)
    p = Patterns(n);

    curr_area1 = p.X_source;
    curr_area2 = p.X_target;
    equals = isequal(p.X_source, p.X_target);

    check = intra_only_overall  &&equals  &&(p.name  ~= "Overall");
    if check
        disp("Skipping; intra_only_overall=" + intra_only_overall + ", equals=" + equals + ", name=" + p.name + ", direction=" + p.directionality);
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
        disp("GOT A")
        result2 = analysis.cca.core(Patterns(end).X_target, ...
            Patterns(end).X_target, 'area1_eq_area2', equals);
        result.b = result2.b;
        result.v = result2.v;
        
    end
    p.cca = result;
    
    % Assign the updated struct back to the original array
    Patterns = nd.set(Patterns, n, p);
    Patterns = nd.set(Patterns, n, p); % TODO: FIX BUG, have to assign twice
    disp(n) 
    assert(~isempty(Patterns(n).cca), "CCA failed for " + Patterns(n).name + " " + Patterns(n).directionality);
end

if remove_inter_from_intra
    disp("Removing inter from intra...")
    % Patterns(1,end).ao = Patterns(1,end).cca.ao;
    % Patterns(1,end).bo = Patterns(1,end).cca.bo;
    disp("GOT B")
    a  = Patterns(1,end).cca.a;
    b  = Patterns(1,end).cca.b;
    br = Patterns(2,end).cca.b;
    ar = Patterns(2,end).cca.a;
    Patterns(1,end).cca.a = munge.remove_via_residual(a, ar);
    Patterns(1,end).cca.b = munge.remove_via_residual(b, br);
    Patterns(1,end).cca.u =  Patterns(1,end).X_source' * Patterns(1,end).cca.a;
    Patterns(1,end).cca.v =  Patterns(2,end).X_target' * Patterns(1,end).cca.b;
end

disp(['CCA analysis completed in ' num2str(toc) ' s']);
