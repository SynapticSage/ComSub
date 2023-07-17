function plotPatternDim(T, genH_name, directionality)
% PLOTPATTERNDIM Plot the dimensionality of the example patterns

    % Create figure
    figName = genH_name + " Example Pattern Dimensionality";
    fig(figName); clf

    % Filter table
    subset = T.genH_name == genH_name & T.directionality == directionality;
    filteredT = T(subset, :);
    assert(~isempty(filteredT))

    % Create gramm plot
    x = categorical(filteredT.patternType);
    y = filteredT.percMax_rrDim;
    g = gramm('x', x, 'y', y, 'color', categorical(filteredT.patternAbstract), 'lightness', categorical(filteredT.control));
    g.geom_jitter('alpha', 0.8, 'width', 0.6);
    g.stat_summary('type', 'sem', 'geom', 'black_errorbar');
    g.set_point_options('base_size', 20);
    g.set_text_options('label_scaling', 1.5, 'base_size', 14);
    g.set_names('x', 'Pattern Type', 'y', '#(Predictive Dims)/#(Max. Pred. Dims)', 'Color', 'Pattern', 'Lightness', 'Treatment/Control');
    g.set_title(genH_name + " Example Pattern Dimensionality" + newline() + directionality);
    g.axe_property('XTickLabelRotation', 35);
    g.draw();

    g.export("file_name", figuredefine("gramm", "patterndim", directionality + "-" + genH_name), "file_type", "svg");
    g.export("file_name", figuredefine("gramm", "patterndim", directionality + "-" + genH_name), "file_type", "pdf");
    g.export("file_name", figuredefine("gramm", "patterndim", directionality + "-" + genH_name), "file_type", "png");

end
