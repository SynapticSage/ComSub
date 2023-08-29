function plotPatternDim(T, genH_name, directionality, varargin)
    % PLOTPATTERNDIM Plot the dimensionality of example patterns
    %

    ip = inputParser();
    ip.addParameter('x', "patternType");
    ip.addParameter('y', "percMax_rrDim");
    ip.addParameter('color', "patternAbstract");
    ip.addParameter('lightness', "control");
    ip.addParameter('pointSize', 8);
    ip.addParameter('displayStatistics', true);
    ip.parse(varargin{:});
    Opt = ip.Results;

    figName = genH_name + " Example Pattern Dimensionality";
    fig(figName); clf

    if ~isempty(directionality)
        subset = T.genH_name == genH_name & T.directionality == directionality;
    else
        subset = T.genH_name == genH_name;
    end
    filteredT = T(subset, :);
    assert(~isempty(filteredT))

    x = categorical(filteredT.patternType);
    y = filteredT.percMax_rrDim;
    g = gramm('x', x, 'y', y, 'color', categorical(filteredT.patternAbstract), 'lightness', categorical(filteredT.control));
    g.geom_jitter('alpha', 0.5, 'width', 0.6, 'height', 0.02);
    g.stat_summary('type', 'sem', 'geom', 'black_errorbar')
    g.set_point_options('base_size', Opt.pointSize)
    g.set_text_options('label_scaling', 1.5, 'base_size', 14);
    g.set_names('x', 'Pattern Type', 'y', '#(Predictive Dims)/#(Max. Pred. Dims)', 'Color', 'Pattern', 'Lightness', 'Treatment/Control');
    g.set_title(upper(genH_name) + " Example Pattern Dimensionality" + newline() + directionality);
    g.axe_property('XTickLabelRotation', 35);
    g.axe_property('YLim', [0, 1]);

    % Conduct statistical tests if displayStatistics is true
    if Opt.displayStatistics
        categories = x;
        uCategories = unique(x);
        nCategories = length(uCategories);
        pValues = [];  % To store p-values
        xPos = [];     % To store x positions for the p-values
        yPos = [];     % To store y positions for the p-values
        xDiff =[];
        pair_name = [];
        for i = progress(1:nCategories-1, 'Title', 'Conducting statistical tests')
            for j = i+1:nCategories
                % Select data for two categories
                a = uCategories(i) == categories;
                b = uCategories(j) == categories;
                data1 = y(a);
                data2 = y(b);

                % Perform Mann-Whitney U test
                [~, p] = ranksum(data1, data2);

                % Store p-value and positions
                pValues = [pValues; p];
                xPos = [xPos; mean([i, j])];
                pair_name = [pair_name; string(uCategories(i)) + "-" + string(uCategories(j))];
                xDiff = [xDiff; i-j];
                if j-i <= 1
                    yPos = [yPos; mean(max([data1, data2]))];
                else
                    yPos = [yPos; -mean(max([data1, data2]))];
                end

                % Plot a line between the two categories if the p-value is
                % significant
                if p < 0.05
                    g.update('x', [i, j], 'y', [max(data1), max(data2)]);
                    g.geom_line();
                    g.set_line_options("styles", "--", "base_size", 2);
                end
            end
        end

        % Create labels for the p-values
        pLabels = compose('%s\np = %.3f', pair_name, pValues);
        g.update('x', xPos, 'y', yPos, 'label', pLabels);
        g.geom_label();
    end

    g.draw();

    set(gcf, 'Position',  get(0, 'Screensize'));
    set(findobj(gcf, 'type', 'line'), 'LineWidth', 2);
    set(findobj(gcf, 'type', 'scatter'), 'AlphaData', 0.332);

    g.export("file_name", figuredefine("gramm", "patterndim", directionality + "-" + genH_name), "file_type", "svg");
    g.export("file_name", figuredefine("gramm", "patterndim", directionality + "-" + genH_name), "file_type", "pdf");
    g.export("file_name", figuredefine("gramm", "patterndim", directionality + "-" + genH_name), "file_type", "png");

    g.draw();

    set(gcf, 'Position',  get(0, 'Screensize'));
    set(findobj(gcf, 'type', 'line'), 'LineWidth', 2);
    set(findobj(gcf, 'type', 'scatter'), 'AlphaData', 0.332);
end

