if ~exist('T', 'var')
    PreambFigs; % Load data
    T(ismember(T.hash, keys), ["genH_name", "preProcess_zscore", "timestamp","animal"])
    % Checkpoint to /Volumes/Ark/commsubspace/
    save(datadefine("tmpT.mat"), 'T')
    T.control = replace(T.control, ["control", "pattern activity"], ["low", "high"]);
end

%                        
% o  /    ,---.,---.,---.
%   /     |---'|---'|---'
%  /      |  \ |  \ |  \ 
% /  o    `   ``   ``   `
%                        

stats = struct();

%% COMPARISON OF FULL MODEL PERF AND PREDICTIVE DIMS
plots.grm.compare_types(T, "power",     "coherence", "hpc-hpc", "axis", 'square');
plots.grm.compare_types(T, "power",     "coherence", "hpc-pfc", "axis", 'square')
plots.grm.compare_types(T, "power",     "wpli",      "hpc-hpc", "axis", 'square');
plots.grm.compare_types(T, "power",     "wpli",      "hpc-pfc", "axis", 'square');
plots.grm.compare_types(T, "coherence", "wpli",      "hpc-hpc", "axis", 'square');
plots.grm.compare_types(T, "coherence", "wpli",      "hpc-pfc", "axis", 'square');

%% SHOW PREDICTIVE DIMENSIONS FOR EACH PATTERN
plots.grm.plotPattern(T, "power",     "hpc-pfc");
plots.grm.plotPattern(T, "power",     "hpc-hpc");
plots.grm.plotPattern(T, "coherence", "hpc-pfc");
plots.grm.plotPattern(T, "coherence", "hpc-hpc");
plots.grm.plotPattern(T, "wpli",      "hpc-pfc");
plots.grm.plotPattern(T, "wpli",      "hpc-hpc");

% ------------------------------------------------------------
% Compare dimensionality for each pattern hpc-hpc vs hpc-pfc
plots.grm.compareField(T, "power",     "field", "rrDim")
plots.grm.compareField(T, "coherence", "field", "rrDim")
plots.grm.compareField(T, "power",     "field", "percMax_rrDim")
plots.grm.compareField(T, "coherence", "field", "percMax_rrDim")
plots.grm.compareField(T, "wpli",      "field", "percMax_rrDim")
plots.grm.compareField(T, "wpli",      "field", "percMax_rrDim")
% plots.grm.compareField(T, "wpli", "field", "rrDim")
% ------------------------------------------------------------
% Test percMax_rrDim
% ------------------------------------------------------------

stats.dim.powcoh.hpchpc = plots.grm.plotWithPermTest(T, "coherence", "power", "hpc-hpc", "field", "percMax_rrdim"); close all
stats.dim.powcoh.hpcpfc = plots.grm.plotWithPermTest(T, "coherence", "power", "hpc-pfc", "field", "percMax_rrdim"); close all

% ------------------------------------------------------------
% Test full_model_performance
% ------------------------------------------------------------
stats.perf.powcoh.hpchpc = plots.grm.plotWithPermTest(T, "coherence", "power", "hpc-hpc", "field", "full_model_performance"); close all
stats.perf.powcoh.hpcpfc = plots.grm.plotWithPermTest(T, "coherence", "power", "hpc-pfc", "field", "full_model_performance"); close all

% ------------------------------------------------------------
% Print stats
% ------------------------------------------------------------

disp("Power vs. coherence hpc-hpc")
struct2table(stats.dim.powcoh.hpcpfc)
disp("Power vs. coherence hpc-pfc")
struct2table(stats.dim.powcoh.hpchpc)
disp("Power vs. coherence hpc-hpc")
struct2table(stats.perf.powcoh.hpcpfc)
disp("Power vs. coherence hpc-pfc")
struct2table(stats.perf.powcoh.hpchpc)

%                                                                       
% ,---.          |                  ,---.          |              o     
% |__. ,---.,---.|--- ,---.,---.    |---|,---.,---.|    ,   .,---..,---.
% |    ,---||    |    |   ||        |   ||   |,---||    |   |`---.|`---.
% `    `---^`---'`---'`---'`        `   '`   '`---^`---'`---|`---'``---'
%                                                       `---'           
%
    %% Cornerhist hpc versus pfc FA qOpt
    clf
    figure(7)
    hpcsubset = T.directionality == 'hpc-hpc';
    pfcsubset = T.directionality == 'hpc-pfc';
    hpcsubset = T(hpcsubset,:);
    pfcsubset = T(pfcsubset,:);
    x = hpcsubset.qOpt;
    y = pfcsubset.qOpt;
    subset =  ~hpcsubset.singularWarning;
    g = gramm(  'subset', subset,...
                'x', x,...
                'y', y,...
                'color',     categorical(hpcsubset.patternAbstract),...
                'lightness', categorical(hpcsubset.control));
    % assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
    g.facet_grid(categorical(hpcsubset.control), categorical(hpcsubset.patternAbstract))
    % g.geom_point('dodge', 0.5, 'alpha', 0.3, 'jitter',0.1);
    g.geom_abline('style','k:');
    g.stat_cornerhist('edges',-4:0.5:4, 'aspect',1.2);
    g.set_point_options('base_size', 10);
    g.set_text_options('label_scaling', 1.5, 'base_size', 14);
    g.set_names('x', "HPC" +newline+ "regional dimensions", ...
        'y', 'PFC regional dimensions', ...
        'Color', 'Pattern', ...
        'Lightness', 'Treatment/Control')
    g.axe_property('XTickLabelRotation',35)
    g.draw()


    %% Cornerhist dimension complexity vs. regression -- PFC more complex but required fewer prediction dimensions
    clf
    figure(8)
    hpcsubset = T.directionality == "pfc-hpc";
    pfcsubset = T.directionality == "pfc-pfc";
    hpcsubset = T(hpcsubset,:);
    pfcsubset = T(pfcsubset,:);
    % x = hpcsubset.qOpt./hpcsubset.rrDim;
    x1 = pfcsubset.qOpt;
    y1 = pfcsubset.rrDim;
    x2 = hpcsubset.qOpt;
    y2 = hpcsubset.rrDim;
    subset = pfcsubset.generateH == "fromSpectral  fromRipTimes" & ~pfcsubset.singularWarning;
    g1 = gramm(  'subset', subset,...
                'x', x1,...
                'y', y1,...
                'color',     categorical(hpcsubset.patternAbstract),...
                'lightness', categorical(hpcsubset.control));


    % assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
    g1.facet_grid(categorical(hpcsubset.control), categorical(hpcsubset.patternAbstract))
    g1.geom_point('dodge', 0.5, 'alpha', 0.3);
    g1.geom_abline('style','k:');
    g1.stat_cornerhist('edges',-4:0.5:4, 'aspect',1.2);
    g1.set_point_options('base_size', 10);
    g1.set_text_options('label_scaling', 1.5, 'base_size', 14);
    g1.set_names('x', 'PFC regional dimensions', ...
        'y', 'PFC predictive dimensions', ...
        'Color', 'Pattern', ...
        'Lightness', 'Treatment/Control')
    g1.axe_property('XTickLabelRotation',35)
    g1.draw()

    %%
    figure(233)
    g2 = gramm(  'subset', subset,...
                'x', x2,...
                'y', y2,...
                'color', categorical(hpcsubset.patternAbstract),...
                'lightness', categorical(hpcsubset.control));
            
    g2.facet_grid(categorical(hpcsubset.control), categorical(hpcsubset.patternAbstract))
    g2.geom_point('dodge', 0.5, 'alpha', 0.3);
    g2.geom_abline('style','k:');
    g2.stat_cornerhist('edges',-4:0.5:4, 'aspect',1.2);
    g2.set_point_options('base_size', 10);
    g2.set_text_options('label_scaling', 1.5, 'base_size', 14);

    g2.set_names('x', 'HPC regional dimensions', ...
        'y', 'HPC predictive dimensions', ...
        'Color', 'Pattern', ...
        'Lightness', 'Treatment/Control')
    g2.axe_property('XTickLabelRotation',35)        
    g2.draw()

% ============================================================
% DIMENSION BAR PLOTS
% ============================================================
% Given Table T
fig('Barplot of percMax_rrDim -- genH == power'); clf;
% Create a gramm object
g = gramm('x', categorical(T.patternAbstractSymbol), 'y', T.percMax_rrDim,...
'color', categorical(T.directionality),...
'subset', T.genH_name == "coherence" & T.control == "pattern activity");
g.set_title('Barplot of percMax_rrDim');

% Set summary statistics for barplot
g.stat_summary('type', 'ci', 'geom', 'bar', 'dodge', 1, 'setylim', true); 

% Aesthetics
g.set_names('x', 'Pattern Abstract Symbol', 'y', 'percMax_rrDim', 'color', 'Directionality');
g.axe_property('XTickLabelRotation', 45); % Rotate x-labels if necessary

% Draw the plot
g.draw();

% Given Table T
fig('animal: Barplot of percMax_rrDim -- genH == power'); clf;
% Create a gramm object
g = gramm('x', categorical(T.patternAbstractSymbol), 'y', T.percMax_rrDim,...
'color', categorical(T.directionality),...
'subset', T.genH_name == "coherence" & T.control == "pattern activity");
g.facet_wrap(categorical(T.animal));
g.set_title('Barplot of percMax_rrDim');

% Set summary statistics for barplot
g.stat_summary('type', 'ci', 'geom', 'bar', 'dodge', 1, 'setylim', true); 

% Aesthetics
g.set_names('x', 'Pattern Abstract Symbol', 'y', 'percMax_rrDim', 'color', 'Directionality');
g.axe_property('XTickLabelRotation', 45); % Rotate x-labels if necessary

% Draw the plot
g.draw();

load RunsSummary
load DetailedRunsSummary
% Convert the timestamp string to datetime format
DetailedRunsSummary.datetime = datetime(DetailedRunsSummary.timestamp, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

% Sort the table based on datetime
DetailedRunsSummary = sortrows(DetailedRunsSummary, 'datetime');

% Get unique directionality values
uniqueDirections = unique(DetailedRunsSummary.directionality);

% Create a color map for the different directionality values
colorMap = lines(length(uniqueDirections));

figure;
hold on;

% Loop over unique directionality values and plot them
for i = 1:length(uniqueDirections)
    % Logical index for current directionality
    idx = DetailedRunsSummary.directionality == uniqueDirections(i);
    
    % Plot the datetime versus percMax_rrDim for this directionality
    markersize = 5;
    if uniqueDirections(i) == "hpc-hpc"
        markersize = 12;
    end
    plot(DetailedRunsSummary.datetime(idx), DetailedRunsSummary.percMax_rrDim(idx), '-o', 'Color', colorMap(i,:), 'MarkerSize', markersize);

end

xlabel('Date');
ylabel('percMax_rrDim');
title('Variation of percMax_rrDim over Time by Directionality');
grid on;

ax = gca;
ax.XTickFormat = 'dd-MMM-yyyy';
xtickangle(45);

% Add a legend
legend(string(uniqueDirections), 'Location', 'best');

% ============================================================


% Get unique animal and directionality values
uniqueAnimals = unique(DetailedRunsSummary.animal);
uniqueDirections = unique(DetailedRunsSummary.directionality);

% Create a color map for the different directionality values
colorMap = lines(length(uniqueDirections));

figure;

% Loop over unique animal values
for a = 1:length(uniqueAnimals)
    % Logical index for current animal
    idx_animal = DetailedRunsSummary.animal == uniqueAnimals(a);
    
    % Create a subplot for this animal
    subplot(length(uniqueAnimals), 1, a);
    hold on;
    
    % Extract the subset of the table for this animal
    subTable = DetailedRunsSummary(idx_animal, :);
    
    % Loop over unique directionality values and plot them
    for i = 1:length(uniqueDirections)
        % Logical index for current directionality
        idx_direction = subTable.directionality == uniqueDirections(i);
        
        markersize = 5;
        if uniqueDirections(i) == "hpc-hpc"
            markersize = 12;
        end
        % Plot the datetime versus percMax_rrDim for this directionality
        plot(subTable.datetime(idx_direction), subTable.percMax_rrDim(idx_direction), '-o', 'Color', colorMap(i,:), 'MarkerSize', markersize);
    end

    xlabel('Date');
    ylabel('percMax_rrDim');
    title(['Animal: ' uniqueAnimals{a}]);
    grid on;

    ax = gca;
    try
    ax.XTickFormat = 'dd-MMM-yyyy';
    catch
    end
    xtickangle(45);
    
    % If it's the last subplot, add a legend
    if a == length(uniqueAnimals)
        legend(string(uniqueDirections), 'Location', 'best');
    end
end

figure;histogram(T.percMax_rrDim(T.directionality=="hpc-hpc" & T.preProcess_zscore==1));hold on;;histogram(T.percMax_rrDim(T.directionality=="hpc-pfc" & T.preProcess_zscore==1))

