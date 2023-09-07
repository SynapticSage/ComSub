function plotIndividSamples(combinedPatternsTable, Option, genH)
%PLOTINDIVIDSAMPLES Plot the performance of individual samples

    combinedPatternsTable = util.table.convertCellstrToString(combinedPatternsTable);

    % Check if genH is provided
    if nargin < 3
        genH = '';
    else
        % Filter the table by genH
        combinedPatternsTable = combinedPatternsTable(combinedPatternsTable.genH == genH, :);
    end

    % Create the plot
    figure('Name','Performance individual samples'); clf
    for i = 1:Option(1).nPatternAndControl
        % col 1 for first 3 and col 2 for last 3
        [row, col] = ind2sub([3,2], i);
        j = sub2ind([2,3], col, row);
        name = Option(1).patternNames(i);
        patternPerformance_pfc = combinedPatternsTable(combinedPatternsTable.direction == "hpc-pfc" & combinedPatternsTable.name == name, :).perf;
        patternPerformance_hpc = combinedPatternsTable(combinedPatternsTable.direction == "hpc-hpc" & combinedPatternsTable.name == name, :).perf;
        subplot(3,2,j)
        p=plot(patternPerformance_pfc);
        set(p,'Color','red','LineWidth',0.25);
        alpha(0.5);
        ylabel("Performance")
        xlabel("Sample")
        hold on
        p=plot(patternPerformance_hpc);
        set(p,'Color','blue','LineWidth',0.25);
        alpha(0.5);
        legend("pfc","hpc")
        title(name)
    end
    % how to interpret nans
    linkaxes(findobj(gcf,'Type','Axes'),'xy')
    ylim(findobj(gcf,'Type','Axes'),[-0.25 0.8])
    if ~isempty(genH)
        sgtitle(genH)
        saveas(gcf,figuredefine("new","examplepredsamps", "exampe_pred_samps_" + genH + ".png"))
        saveas(gcf,figuredefine("new","examplepredsamps", "exampe_pred_samps_" + genH + ".pdf"))
    else
        saveas(gcf,figuredefine("new","examplepredsamps", "exampe_pred_samps.png"))
        saveas(gcf,figuredefine("new","examplepredsamps", "exampe_pred_samps.pdf"))
    end
end
