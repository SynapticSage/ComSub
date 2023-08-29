function barplot_property(T, y_var_name, varargin)
    % Input parsing
    p = inputParser;
    addRequired(p, 'T');
    addRequired(p, 'y_var_name', @ischar);
    addParameter(p, 'ylim_percentiles', [0.05, 0.95], @isnumeric);
    addParameter(p, 'individual_ylims_for_animals', false, @islogical);
    addParameter(p, 'facet_by_animal', false, @islogical);
    
    parse(p, T, y_var_name, varargin{:});
    
    % Extract the parsed inputs
    T = p.Results.T;
    y_var_name = p.Results.y_var_name;
    ylim_percentiles = p.Results.ylim_percentiles;
    individual_ylims_for_animals = p.Results.individual_ylims_for_animals;
    facet_by_animal = p.Results.facet_by_animal;

    % The rest of the function remains unchanged...

    clf;

    % Extract y data
    y_data = T.(y_var_name);

    if individual_ylims_for_animals
        uniqueAnimals = unique(T.animal);
    else
        y_lower_quantile = quantile(y_data, ylim_percentiles(1));
        y_upper_quantile = quantile(y_data, ylim_percentiles(2));
    end

    % Create a gramm object
    g = gramm('x', categorical(T.patternAbstractSymbol), 'y', y_data, 'color', categorical(T.directionality));
    g.set_title(['Barplot of ', y_var_name]);
    g.stat_summary('type', 'ci', 'geom', {'bar', 'errorbar'}, 'dodge', 1, 'setylim', true);
    g.set_names('x', 'Pattern Abstract Symbol', 'y', y_var_name, 'color', 'Directionality');
    g.axe_property('XTickLabelRotation', 45); % Rotate x-labels if necessary

    % Faceting
    if facet_by_animal
        g.facet_wrap({categorical(T.animal), categorical(T.genH_name)});
    else
        g.facet_wrap(categorical(T.genH_name));
    end

    % Set global y-lims if not setting individually
    if ~individual_ylims_for_animals
        g.axe_property('YLim', [y_lower_quantile, y_upper_quantile]);
    end

    g.set_text_options('Interpreter', 'tex', 'base_size', 14); % Set text options

    % Draw the plot
    g.draw();

    % Loop through animals to set individual y-lims if flag is set
    if individual_ylims_for_animals
        for a = 1:length(uniqueAnimals)
            subplotHandle = subplot(g.facet_axes_handles(a));
            animalData = y_data(T.animal == uniqueAnimals(a));
            y_lower_quantile = quantile(animalData, ylim_percentiles(1));
            y_upper_quantile = quantile(animalData, ylim_percentiles(2));
            ylim(subplotHandle, [y_lower_quantile, y_upper_quantile]);
        end
    end

    % Exporting the figure
    if ~exist(figuredefine("gramm"), 'dir')
        mkdir(figuredefine("gramm"))
    end
    print(gcf, figuredefine("gramm", [y_var_name, "_facet"]), '-dpdf', '-vector');
end

