function [h_corrdiff,p_corrdiff] = plotCfByDirection(FigDat, Patterns, varargin)

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
ip.addParameter('normalization', 'count', @(x) ischar(x) || isstring(x));
ip.addParameter('example', 1, @(x) isnumeric(x));
ip.parse(varargin{:});
Opt = ip.Results;

% How about just put all the patterns together?
fig 'Overall directional correlation'
subplot(2,1,1)
ax1 = nexttile;
hist_withpfc = histogram(all_pairs_withhpc, 'Normalization',Opt.normalization);
ylabel("Pairs")
title ("HPC-HPC")
example_mean_corr_withhpc = FigDat.mean_corrwithhpc(Opt.example);
example_std_corr_withhpc  = FigDat.std_corrwithhpc(Opt.example);

hold on
x=[example_mean_corr_withhpc, example_mean_corr_withhpc];
y=[0 max(hist_withpfc.Values)];
lineObject=line(x,y);
lineObject.LineStyle = ':'; % Make line dotted
lineObject.LineWidth = 2;  % Thicken the line
lineObject.Color = 'black'; % Color it black
subplot(2,1,2)
ax2 = nexttile;
hist_hp = histogram(all_pairs_withpfc,'Normalization',Opt.normalization);
ylabel("Pairs")
title ("HPC-PFC")
example_mean_corr_withpfc = FigDat.mean_corrwithpfc(Opt.example);
example_std_corr_withpfc  = FigDat.std_corrwithpfc(Opt.example);

x=[example_mean_corr_withpfc,example_mean_corr_withpfc];
y=[0 max(hist_hp.Values)]
lineObject=line(x,y);
lineObject.LineStyle = ':'; % Make line dotted
lineObject.LineWidth = 2;  % Thicken the line
lineObject.Color = 'black'; % Color it black
xlabel("Pairwise correlation")
linkaxes([ax1,ax2],'x');
[h_corrdiff,p_corrdiff] = kstest2(all_pairs_withpfc,all_pairs_withhpc);
%%
% print stats
formatSpec1 = "%s: %0.3f±%0.3f";

sprintf(formatSpec1,Patterns(1,1).directionality,mean(mean_corrwithhpc),mean(std_corrwithhpc))
sprintf(formatSpec1,Patterns(2,1).directionality,mean(mean_corrwithpfc),mean(std_corrwithpfc))
disp(p_corrdiff)

examp_str = sprintf("Example %d, %s ", Opt.example, Option(1).patternNamesFull(Opt.example));
fullfigname = fullfile(codedefine,'figures', 'new','cofire',...
    examp_str + string(Opt.normalization) + "cofiring per pattern" + Opt.figAppend);
saveas(gcf, fullfigname + ".png")
saveas(gcf, fullfigname + ".pdf")
savefig(gcf, fullfigname + ".fig")

