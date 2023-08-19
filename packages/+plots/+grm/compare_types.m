function compare_types(T, genH_name1, genH_name2, directionality, varargin)
% PLOTCOMPARISON Plot comparison of two genH values
%   PLOTCOMPARISON(T, genH_name1, genH_name2) plots the comparison of two
%   genH values in table T. genH_name1 and genH_name2 must be one of the
%   values in the genH column of T. 

    ip = inputParser();
    ip.addParameter('point_size', 3);
    ip.addParameter('investigate', ["percMax_rrDim", "full_model_performance"]);
    ip.addParameter('investNames', ["Predictive Dims", "Model Performance"]);
    ip.addParameter('axis',''); % if square, make the axes square
    ip.parse(varargin{:});
    point_size = ip.Results.point_size;
    investigate = ip.Results.investigate;
    investNames = ip.Results.investNames;
    axshape = ip.Results.axis;

    % Create figure
    figName = genH_name1 + " vs " + genH_name2 + " <" + directionality + ">";
    fig(figName); clf

    % Filter table
    subset1 = T.genH_name == genH_name1 & T.directionality == directionality;
    subset2 = T.genH_name == genH_name2 & T.directionality == directionality;
    filteredsubset = T(subset1, :);
    coherencesubset = T(subset2, :);
    assert(~isempty(filteredsubset))
    assert(~isempty(coherencesubset))
    assert(height(filteredsubset) == height(coherencesubset))

    % Create gramm plots
    for i = 1:length(investigate)
        varName = investigate(i);
        xlabelSuffix = investNames(i);
        M = max([abs(filteredsubset.(varName)); abs(coherencesubset.(varName))]);
        corner_kws = {'edges', -M:(1/10 * M):M, 'aspect', 0.6, 'location', M/2, 'fill', 'transparent', 'normalization', 'countdensity'};
        g(1,i) = createPlot(filteredsubset, coherencesubset, varName, corner_kws, point_size, xlabelSuffix, genH_name1, genH_name2);
    end

    if strcmp(axshape, 'square')
        g.axe_property('DataAspectRatio', [1 1 1]);
    end

    % Draw the plots
    g.draw();
    % for i = 1:numel(g)
    %     for ax = g(i).facet_axes_handles
    %         axis(ax, 'square');
    %     end
    % end
    set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure.

    % Save the figure
    g.export("file_name", figuredefine("gramm", "compare_perf_and_dim", figName + ".png"));
    g.export("file_name", figuredefine("gramm", "compare_perf_and_dim", figName + ".pdf"));

end

function g = createPlot(filteredsubset, coherencesubset, varName, corner_kws, point_size, xlabelSuffix, genH_name1, genH_name2)
    filteredsubset.control = replace(filteredsubset.control, ["control", "pattern activity"], ["low", "high"]);
    coherencesubset.control = replace(coherencesubset.control, ["control", "pattern activity"], ["low", "high"]);
    x = filteredsubset.(varName);
    y = coherencesubset.(varName);
    g = gramm('x', x, 'y', y, 'color', categorical(filteredsubset.patternAbstract), 'lightness', categorical(filteredsubset.control));
    g.facet_grid(categorical(filteredsubset.control), categorical(filteredsubset.patternAbstract));
    g.geom_point('dodge', 0.5, 'alpha', 0.3);
    g.geom_abline('style', 'k:');
    g.stat_cornerhist(corner_kws{:});
    g.set_point_options('base_size', point_size);
    g.set_text_options('label_scaling', 1.2, 'base_size', 14);
    g.set_names('x', genH_name1 + newline + xlabelSuffix, 'y', genH_name2 + newline + xlabelSuffix, 'Color', 'Pattern', 'Lightness', 'Treatment/Control');
end

