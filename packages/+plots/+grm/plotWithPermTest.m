function stats = plotWithPermTest(T, genH_name1, genH_name2, directionality, varargin)
    % PLOTCORNERHIST Plot the corner histogram of the example patterns
    %
    %   PLOTCORNERHIST(T, genH_name1, genH_name2, directionality) plots the
    %   corner histogram of the example patterns in table T. genH_name1 and
    %   genH_name2 must be one of the values in the genH column of T. The 
    %   directionality must be either "hpc-hpc" or "hpc-pfc".

    ip = inputParser();
    ip.addParameter('pointSize', 10);
    ip.addParameter('alpha', 0.1);
    ip.addParameter('width', 0.35);
    ip.addParameter('height', 0.35);
    ip.addParameter('numPermutations', 1000);
    ip.addParameter('showPermutationTest', true);
    ip.addParameter('field', "percMax_rrDim", @isstring);
    ip.parse(varargin{:});
    Opt = ip.Results;

    assert(genH_name1 ~= genH_name2, "genH_name1 and genH_name2 must be different");

    % Filter table
    subset1 = T.genH_name == genH_name1 & T.directionality == directionality;
    subset2 = T.genH_name == genH_name2 & T.directionality == directionality;
    filteredsubset1 = T(subset1, :);
    filteredsubset2 = T(subset2, :);
    assert(~isempty(filteredsubset1))
    assert(~isempty(filteredsubset2))
    disp("Number of example patterns in " + genH_name1 + ": " + num2str(height(filteredsubset1)));
    disp("Number of example patterns in " + genH_name2 + ": " + num2str(height(filteredsubset2)));

    % Set corner histogram parameters
    M = max([abs(filteredsubset1.(Opt.field)); ... 
             abs(filteredsubset2.(Opt.field))]);
    m = min([abs(filteredsubset1.(Opt.field)); ... 
             abs(filteredsubset2.(Opt.field))]);
    corner_kws = {'edges', -M:(1/20 * M):M, 'aspect', 1, 'location', M/2, 'fill', 'transparent', 'normalization', 'countdensity'};

    % Create figure
    f = fig("permtest: " + genH_name1 + " vs " + genH_name2 + " (" + directionality + ")"); clf;
    
    % Create gramm object for control group
    g = gramm('x', filteredsubset1.(Opt.field), 'y', filteredsubset2.(Opt.field), 'subset', filteredsubset1.control == "control");
    g.geom_jitter('alpha', Opt.alpha, 'width', Opt.width, 'height', Opt.height);
    g.stat_cornerhist(corner_kws{:});
    g.set_point_options('base_size', Opt.pointSize);
    g.set_text_options('label_scaling', 1.5, 'base_size', 10);
    g.set_names('x', genH_name1 + " Predictive Dims",...
    'y', genH_name2 + " Predictive Dims",...
    'Color', 'Pattern', 'row', '', 'Lightness', 'Treatment/Control');
    g.set_color_options('chroma', 0);
    g.set_text_options('interpreter', 'latex', 'base_size', 10);
    g.axe_property('XLim', [-m M], 'YLim', [-m M]);
    g.draw();

    % Update gramm object for pattern activity group
    g.update('subset', filteredsubset1.control ~= "control", 'color', categorical(filteredsubset1.patternAbstract));
    g.set_color_options();
    g.stat_cornerhist(corner_kws{:});
    g.geom_abline('style', 'k:');
    g.geom_jitter('alpha', 0.10, 'width', Opt.width, 'height', Opt.height);
    g.draw();

    p           = nan(size(unique(filteredsubset1.patternAbstract)));
    obsdiffs    = nan(size(unique(filteredsubset1.patternAbstract)));
    effectsizes = nan(size(unique(filteredsubset1.patternAbstract)));
    patAbs      = string(nan(size(unique(filteredsubset1.patternAbstract))));

    % Conduct permutation tests if showPermutationTest is true
    if Opt.showPermutationTest
        [G, patAb] = findgroups(filteredsubset1.patternAbstract);
        for gg = unique(G)'
            a = filteredsubset1(G == gg & ...
            filteredsubset1.control == "pattern activity", :);
            b = filteredsubset2(G == gg & ...
            filteredsubset2.control == "pattern activity",:);
            patA = a.patternAbstract(1);
            a = a.(Opt.field);
            b = b.(Opt.field);
            % a = a(randperm(length(a)));
            % b = b(randperm(length(b)));
            t_pat = a - b;
            a = filteredsubset1(G == gg & filteredsubset1.control == "control", :).(Opt.field);
            b = filteredsubset2(G == gg & filteredsubset2.control == "control",:).(Opt.field);
            % a = a(randperm(length(a)));
            % b = b(randperm(length(b)));
            t_ctrl = filteredsubset1(G == gg & filteredsubset1.control == "control", :).(Opt.field) - filteredsubset2(G==gg & filteredsubset2.control=="control",:).(Opt.field);
            ploton = 0;
            [p(gg), obsdiffs(gg), effectsizes(gg)] = ... 
            permutationTest(t_pat, t_ctrl, Opt.numPermutations, 'plotresult', ploton, 'meanfunc', @nanmedian);
            patAbs(gg) = patA;
            if ploton
                title(patAb(gg));
                set(gca, 'YScale', 'log');
            end
            
        end
    end

    % Add median to each corner_hist
    [G, patAb, ctrl] = findgroups(filteredsubset1.patternAbstract, filteredsubset1.control);
    [~,~,uPatAb] = unique(patAb);
    [~,~,uCtrl] = unique(ctrl);
    med = splitapply(@median, filteredsubset1.(Opt.field) - filteredsubset2.(Opt.field), G); % get median of each group
    for gg = unique(G)'
        corner_hist_handle = g.results.stat_cornerhist(uPatAb(gg)).child_axe_handle;
        if ctrl(gg) == "control"
            corner_hist_color = [0.5, 0.5, 0.5];
        else
            corner_hist_color = g.results.stat_cornerhist(uPatAb(gg)).bar_handle.FaceColor / 2;
        end
        Y = ylim(corner_hist_handle);
        line(corner_hist_handle, [med(gg) med(gg)], Y*1.1, 'Color', corner_hist_color, 'LineWidth', 3, 'LineStyle', ':');
    end

    stats.p = p;
    stats.obsdiffs = obsdiffs;
    stats.effectsizes = effectsizes;
    stats.patAbs = patAbs;
    stats.genH_name1 = repmat(genH_name1, size(p));
    stats.genH_name2 = repmat(genH_name2, size(p));
    stats.field = repmat(Opt.field, size(p));

end

