function plotDimensionalRemoval_perPatternbyDirection(rt, varargin)
% plotDimensionalRemoval_perPatternbyDirection(rt, varargin)

ip = inputParser;
ip.addParameter('figAppend', []);
ip.parse(varargin{:});
Opt = ip.Results;
Opt.figAppend = string(Opt.figAppend);
eb = 'sem';

eb = @(y)([trimmean(y,2.5);bootci(500,{@(ty)trimmean(ty,2.5),y},'alpha',0.05)]);
non_high = ["control", "mid"];
not_control = ~contains(rt.removePattern, non_high) & ~contains(rt.basePattern,   non_high);
control = ~contains(rt.removePattern, "control") & contains(rt.basePattern,   "control");

figFolder = fullfile(figuredefine, 'dimensionRemoval');
if ~exist(figFolder, "dir")
    mkdir(figFolder);
end
skws = {};

[groups, genHs] = findgroups(rt.method);
uGroups = unique(groups);
for G = uGroups'

    RT = rt(groups == G,:);
    figc("dimensionRemoval perPattern" + genHs(G));
    g = gramm(...
        'x',RT.dimensionRemoved, ...
        'marker', categorical(RT.targetArea),...
        'y', RT.performance, ...
        'color', categorical(RT.sameDirectionLabel),'subset',~RT.sameDirection);
    g = g.facet_grid([],categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar', 'type', eb);
    g = g.stat_summary('geom','area', 'type', eb);
    g = g.set_color_options('chroma',0,'lightness',30);
    g.set_text_options("interpreter",'latex');
    g = g.set_names('x','dims removed',...
        'y','peformance',...
        'column','Interaction',...
        'row','Pattern', ...
        'color','Brain area removed');
    snapnow;
    g=g.update('subset',RT.sameDirection);
    %g = g.facet_grid([],categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar', 'type', eb);
    g = g.stat_summary('geom','area', 'type', eb);
    g = g.set_color_options(); % Restore default color
    g.draw();
    warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
    sgtitle("dimensionRemoval perPattern" + genHs(G));
    g.export('file_name', figuredefine("dimensionRemoval", get(gcf,'Name') + ".svg"),skws{:})
    g.export('file_name', figuredefine("dimensionRemoval", get(gcf,'Name') + ".png"),skws{:})
end

tmp = contains(rt.basePattern,'control') & ~contains(rt.basePattern,"mid");
rt.control = repmat("",size(rt,1),1);
rt.control(tmp) = "low";
tmp = ~contains(rt.basePattern,'control') & ~contains(rt.basePattern,"mid");
rt.control(tmp) = "high";
tmp = contains(rt.basePattern,"mid");
rt.control(tmp) = "mid";

[groups, genHs] = findgroups(rt.method);
uGroups = unique(groups);
for G = uGroups'

    RT = rt(groups == G,:);
    figc("dimensionRemoval perPattern" + genHs(G));
    g = gramm(...
        'x',RT.dimensionRemoved, ...
        'marker', categorical(RT.targetArea),...
        'y', RT.performance, ...
        'color', categorical(RT.sameDirectionLabel),'subset',~RT.sameDirection);
    g = g.facet_grid(categorical(RT.control),categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar', 'type', eb);
    g = g.stat_summary('geom','area', 'type', eb);
    g = g.set_color_options('chroma',0,'lightness',30);
    g.set_text_options("interpreter",'latex');
    g = g.set_names('x','dims removed',...
        'y','peformance',...
        'column','Interaction',...
        'row','Pattern', ...
        'color','Brain area removed');
    snapnow;
    g=g.update('subset',RT.sameDirection);
    %g = g.facet_grid([],categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar', 'type', eb);
    g = g.stat_summary('geom','area', 'type', eb);
    g = g.set_color_options(); % Restore default color
    % set xlims
    g.axe_property('XLim',[0 10]);
    g.draw();
    warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
    sgtitle("dimensionRemoval perPattern" + genHs(G) + newline + Opt.figAppend + newline + newline);
    export_name = replace(get(gcf,'Name'),newline, "_");
    g.export('file_name', figuredefine("dimensionRemoval", export_name), skws{:})
    g.export('file_name', figuredefine("dimensionRemoval", export_name), skws{:})
end
