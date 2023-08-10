function savedat(Patterns_overall, Spk, behavior, Option, figAppend)
    spikeRateMatrix = Spk.spikeRateMatrix;
    areaPerNeuron   = Spk.areaPerNeuron;
    srm_time        = Spk.timeBinMidPoints;
    [U, V] = deal(Patterns_overall(end).cca.u, Patterns_overall(end).cca.v);
    uv_time = Patterns_overall(end).X_time;
    % Convert behavior table to pandas DataFrame
    behavior_dict = table2struct(behavior, 'ToScalar', true);
    behavior_names = behavior.Properties.VariableNames;
    if ~exist(figuredefine("data"), 'dir')
        mkdir(figuredefine("data"));
    end
    file = figuredefine("data", figAppend + ".mat");
    disp("Saving data to " + file);
    save(file, 'spikeRateMatrix', 'areaPerNeuron', 'srm_time', ...
    'U', 'V', 'uv_time', 'Option', 'behavior_dict', 'behavior_names');
end

