function regionalPerf(combinedPatternsTable, singlePredTable, varargin)
% regionalPerf

const = option.constants();
ip = inputParser();
ip.addParameter('select_genH', [], @(x) iscellstr(x) || isstring(x) || ischar(x));
ip.parse(varargin{:});
Opt = ip.Results;

combinedPatternsTable.genH = string(combinedPatternsTable.genH);
Opt.select_genH = string(Opt.select_genH);
if ~isempty(Opt.select_genH)
    combinedPatternsTable = combinedPatternsTable(ismember(combinedPatternsTable.genH, Opt.select_genH), :);
    assert(~isempty(combinedPatternsTable), 'combinedPatternsTable is empty');
end

r_square_withhpc = combinedPatternsTable.perf(combinedPatternsTable.direction == "hpc-hpc");
r_square_withpfc = combinedPatternsTable.perf(combinedPatternsTable.direction == "hpc-pfc");
assert(~isempty(r_square_withhpc), 'r_square_withhpc is empty');
assert(~isempty(r_square_withpfc), 'r_square_withpfc is empty');
if ismember('Var8', singlePredTable.Properties.VariableNames) && ~ismember('perf', singlePredTable.Properties.VariableNames)
    singlePredTable.Properties.VariableNames{'Var8'} = 'perf';
end
mean_withhpc = mean(r_square_withhpc(intersect(~isinf(r_square_withhpc), ~isnan(r_square_withhpc)))); % mean of r2 with hpc
std_withhpc  = std(r_square_withhpc(intersect(~isinf(r_square_withhpc),  ~isnan(r_square_withhpc)))); % std of r2 with hpc
mean_withpfc = mean(r_square_withpfc(intersect(~isinf(r_square_withpfc), ~isnan(r_square_withpfc)))); % mean of r2 with pfc
std_withpfc  = std(r_square_withpfc(intersect(~isinf(r_square_withpfc),  ~isnan(r_square_withpfc)))); % std of r2 with pfc
if ~isempty(singlePredTable)
    singlePredTable.genH = string(singlePredTable.genH);
    singlePredTable = singlePredTable(ismember(singlePredTable.genH, Opt.select_genH), :);
    single_pred_with_hpc = singlePredTable.perf(singlePredTable.direction == "hpc-hpc");
    single_pred_with_pfc = singlePredTable.perf(singlePredTable.direction == "hpc-pfc");
    % Compute medians
    median_singlehh = median(single_pred_with_hpc(~isnan(single_pred_with_hpc))); % median of single neuron prediction with hpc
    median_singlehp = median(single_pred_with_pfc(~isnan(single_pred_with_pfc))); % median of single neuron prediction with pfc
end

%%
[m, M] = deal(-0.1, 0.4);
ax1=subplot(2,1,1);
% constrain to -1 to 1
r_square_withhpc(r_square_withhpc > M) = M;
r_square_withhpc(r_square_withhpc < m) = m;
hist_withhpc = histogram(r_square_withhpc,50, 'FaceColor', const.hpccolor);
ylabel("Data Sets")
title ("source predicting HPC targets")
lineObject=line([mean_withhpc,mean_withhpc],[0 max(hist_withhpc.Values)]);
lineObject.LineStyle = ':'; % Make line dotted
lineObject.LineWidth = 2;  % Thicken the line
lineObject.Color = 'black'; % Color it black
if ~isempty(singlePredTable)
    lineObject2 = line([median_singlehh,median_singlehh],[0 max(hist_withhpc.Values)]);
    lineObject2.LineWidth = 3;
    lineObject2.Color = 'blue';
end
%%
ax2=subplot(2,1,2)
% Constrain to -1 to 1
r_square_withpfc(r_square_withpfc > M) = M;
r_square_withpfc(r_square_withpfc < m) = m;
hist_hp = histogram(r_square_withpfc,50, 'FaceColor', const.pfccolor); % histogram of 100 bins
% prediction performance from hpc to pfc cells
ylabel("Data Sets")
title ("source predicting PFC targets")
xlabel("Performance")
linkaxes([ax1,ax2],'x');
lineObject=line([mean_withpfc,mean_withpfc],[0 max(hist_hp.Values)]);
lineObject.LineStyle = ':'; % Make line dotted
lineObject.LineWidth = 2;  % Thicken the line
lineObject.Color = 'black'; % Color it black
if ~isempty(singlePredTable)
    lineObject2 = line([median_singlehp,median_singlehp],[0 max(hist_hp.Values)]);
    lineObject2.LineWidth = 3;  % Thicken the line
    lineObject2.Color = 'blue'; % Color it black
end
