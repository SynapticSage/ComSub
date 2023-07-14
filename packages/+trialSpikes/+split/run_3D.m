function out = run(Spk, Option)
% trialSpikes.split.run
%
% Input:
%   Spk: trialSpikes object
%   Option: struct
%       .binsToMatchFR: number of bins to match FR
%       .sourceArea: 'hpc' or 'pfc'
%
% Output:
%   out: struct
%       .source: source spikes
%       .target: target spikes
%       .nSource: number of source cells
%       .nTarget: number of target cells
%       .source_index: index of source cells
%       .target_index: index of target cells

    out.sourceArea = Option.sourceArea;

    if strcmp(Option.sourceArea, "CA1")
        source = Spk.hpc.T;
        other_target = Spk.pfc.T;
        out.directionality = ["hpc-hpc", "hpc-pfc"];
    elseif strcmp(Option.sourceArea, "PFC")
        source = Spk.pfc.T;
        other_target = Spk.hpc.T;
        out.directionality = ["pfc-pfc", "pfc-hpc"];
    else
        error("Invalid source area");
    end

    [source, target, nSource, nTarget, source_index, target_index] ...
        = trialSpikes.split.twoWayFRMatch(Spk.hpc.T, Spk.pfc.T, ...
                      Spk.hpc.FR, Spk.pfc.FR, Option.binsToMatchFR);
    out.source = source;
    out.target = target;
    out.nSource = nSource;
    out.nTarget = nTarget;
    out.source_index = source_index;
    out.target_index = target_index;
    out.source_area = Option.sourceArea;

end