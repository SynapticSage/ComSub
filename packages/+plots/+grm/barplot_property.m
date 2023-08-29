function barplot_property(T, y_var_name, varargin)
% barplot_property(T, y_var_name, varargin)
%
% Plots a barplot of the given y_var_name from the table T.
%
% Inputs:
%   T - A table containing the data to plot
%   y_var_name - The name of the variable to plot on the y-axis
%
% Optional Inputs:
%   ylim_percentiles - A 1x2 vector containing the lower and upper
%       percentiles to use for the y-axis limits. Default is [0.05, 0.95].
%   individual_ylims_for_animals - A boolean indicating whether to set
%       individual y-limits for each animal. Default is false.
%   facet_by_animal - A boolean indicating whether to facet the plot by
%       animal. Default is false.

    % Input parsing
    p = inputParser;
    addParameter(p, 'ylim_percentiles', [0.05, 0.95], @isnumeric);
    addParameter(p, 'individual_ylims_for_animals', false, @islogical);
    addParameter(p, 'facet_by_animal', false, @islogical);
    addParameter(p, 'switch_facet', false, @islogical);
    
    p.parse(varargin{:});
    
    % Extract the parsed inputs
    ylim_percentiles = p.Results.ylim_percentiles;
    individual_ylims_for_animals = p.Results.individual_ylims_for_animals;
    facet_by_animal = p.Results.facet_by_animal;
    switch_facet = p.Results.switch_facet;

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
        if switch_facet
            g.facet_grid(categorical(T.genH_name), categorical(T.animal));
        else
            g.facet_grid(categorical(T.animal), categorical(T.genH_name));
        end
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
    if facet_by_animal
        facet = "_facet_by_animal_and_genH";
    else
        facet = "_facet_by_genH";
    end
    if ~exist(figuredefine("gramm"), 'dir')
        mkdir(figuredefine("gramm"))
    end
    print(gcf, figuredefine("gramm", "barplot_" + y_var_name + "_facet" + facet), '-dpdf', '-vector','-fillpage');
end
