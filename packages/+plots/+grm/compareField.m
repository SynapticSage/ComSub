function compareField(T, genH_name, varargin)
%PLOTCORNERHIST Plot a corner histogram of the predictive dimensions of

    ip = inputParser();
    ip.addParameter('field', "maxPerc_rrDim", @isstring);
    ip.addParameter('axis', "square", @isstring);
    ip.addParameter('alpha', 0.25, @isnumeric);
    ip.addParameter('base_size', 5, @isnumeric);
    ip.parse(varargin{:});
    compare = ip.Results.field;
    axshape = ip.Results.axis;
    Opt = ip.Results;

    % Define subsets
    hpcsubset = T(T.directionality == "hpc-hpc" & T.genH_name == genH_name,:);
    pfcsubset = T(T.directionality == "hpc-pfc" & T.genH_name == genH_name,:);

    % Define variables for plotting
    x = hpcsubset.(compare);
    y = pfcsubset.(compare);


    g=runplot(x, y, hpcsubset, Opt, genH_name);

    % Save the plot
    if ~exist(figuredefine("gramm", compare))
        mkdir(figuredefine("gramm", compare));
    end
    g.export("file_name", figuredefine("gramm", compare, "compare-dir-" + genH_name + "_cornerhist"), "file_type", "pdf");
    g.export("file_name", figuredefine("gramm", compare, "compare-dir-" + genH_name + "_cornerhist"), "file_type", "png");  

    % Draw the plot
    try
        runplot(x, y, hpcsubset, Opt, genH_name);
    catch ME
    end


function g= runplot(x, y, hpcsubset, Opt, genH_name)

    compare = Opt.field;
    axshape = Opt.axis;
    % Calculate the maximum absolute value of x and y
    M = max([abs(x); abs(y)]);
    m = min([abs(x); abs(y)]);
    j = (M-m)*0.01;

    % Set up figure
    clf
    f = fig("compareField: " + compare + " " + genH_name);
    corner_kws = {'edges', -M:(1/20 * M):M, 'aspect', 1, 'location', [M/2], 'fill', 'transparent', 'normalization', 'countdensity'};

    % Create gramm plot
    g = gramm('x', x, 'y', y, 'subset', hpcsubset.control == "low");
    g.facet_grid(categorical(hpcsubset.patternAbstractSymbol), []);
    g.geom_jitter('alpha', Opt.alpha, 'width', j, 'height', j);
    g.stat_cornerhist(corner_kws{:});
    g.set_point_options('base_size', Opt.base_size * 1.25);
    g.set_text_options('label_scaling', 1.5, 'base_size', 10);
    g.set_names('x', "HPC" + newline + compare, 'y', "PFC" + newline + compare, 'Color', 'Pattern', 'row', '', 'Lightness', 'Treatment/Control');
    g.set_color_options('chroma', 0);
    g.set_text_options('interpreter', 'latex', 'base_size', 10);

    % Update gramm plot
    g.update('subset', hpcsubset.control == "pattern activity", 'color', categorical(hpcsubset.patternAbstract));
    % set to gray
    g.set_point_options('base_size', Opt.base_size);
    g.set_color_options();
    g.stat_cornerhist(corner_kws{:});
    g.geom_abline('style', 'k:');
    g.geom_jitter('alpha', Opt.alpha, 'width', j, 'height', j);
    if axshape == "square"
        g.axe_property('DataAspectRatio', [1 1 1]);
    end

    % Draw the plot
    g.draw();
    set(gcf, 'Position', get(0, 'ScreenSize'));

