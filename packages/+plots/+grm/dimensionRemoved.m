function dimensionRemoved(rt, genH, varargin)
% DIMENSIONREMOVED Plot the performance of the model when removing

ip = inputParser();
ip.addParameter('figAppend', "", @(x) ischar(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;
Opt.figAppend = string(Opt.figAppend);
% errorbar method
% eb = "sem";
eb = @(y)([trimmean(y,2.5);bootci(500,{@(ty)trimmean(ty,2.5),y},'alpha',0.05)]);
control_pattern = ["control", "mid"];

if ~exist(figuredefine("dimensionRemoval"))
    mkdir(figuredefine("dimensionRemoval"))
end

rtnew  = rt(rt.method == genH,:);
assert(~isempty(rtnew), "No data for method " + genH);
not_control = ~contains(rtnew.removePattern, control_pattern) & ~contains(rtnew.basePattern,   control_pattern);

figure;
g = gramm(...
    'x', rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color',     categorical(rtnew.removePattern), ...
    'lightness', categorical(rtnew.sameDirectionLabel),...
    'linestyle', categorical(rtnew.sameDirectionLabel),...
    'subset', not_control);
g = g.facet_grid(categorical(rtnew.basePatternLabel), ...
                 categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5,'type',eb);
g = g.stat_summary('geom','area','type',eb);
g.set_text_options("interpreter",'latex');
g = g.set_names('x','dims removed',...
                'y','peformance',...
                'column','Interaction',...
                'row','Pattern', ...
                'color','Pattern dims removed',...
                'lightness',"Remove Same/Different"+ newline+"Pred. Target",...
                'linestyle',"Remove Same/Different"+ newline+"Pred. Target");
g.axe_property('XLim',[0 6]);
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval",genH + "-remove dims from same+diff - zoom" + Opt.figAppend),'file_type',"svg");
g.export('file_name',figuredefine("dimensionRemoval",genH + "-remove dims from same+diff" + Opt.figAppend),'file_type',"pdf");
poststeps(g)

%%  REMOVE PATTERN FROM SAME TARGET AREA ONLY
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',rtnew.sameDirection & not_control);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5,'type',eb);
% g = g.stat_summary('geom','area','type',eb);
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of similar target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval",genH + "-remove dims from same target area only" + Opt.figAppend),'file_type',"pdf");
g.export('file_name',figuredefine("dimensionRemoval",genH + "-remove dims from same target area only - zoom" + Opt.figAppend),'file_type',"svg");
poststeps(g)

% REMOVE PATTERN FROM SAME TARGET AREA ONLY, zoom into [0,6]
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',rtnew.sameDirection & not_control);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5,'type','sem','type',eb);
% g = g.stat_summary('geom','area','type',eb);
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of similar target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
set(g.facet_axes_handles,'xlim',[0 6])
steps()

g.export('file_name',figuredefine("dimensionRemoval", genH + "-remove dims from same target area only - zoom" + Opt.figAppend),'file_type',"svg");
g.export('file_name',figuredefine("dimensionRemoval", genH + "-remove dims from same target area only" + Opt.figAppend),'file_type',"pdf");
poststeps(g)

%% REMOVE PATTERN FROM DIFFERENT TARGET AREA
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',~rtnew.sameDirection & not_control);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','type',eb);
g = g.stat_summary('geom','area','type',eb);
g = g.set_color_options('lightness',200);
g.set_text_options("interpreter",'latex');
g = g.set_names('x', 'dims removed',...
                'y', 'peformance',...
                'column', 'Interaction',...
                'row', 'Pattern', ...
                'color', 'Pattern dims removed');
g=g.set_title('Removing different only');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval", genH + "-remove dims from different target area only" + Opt.figAppend),       'file_type',"pdf");
g.export('file_name',figuredefine("dimensionRemoval", genH + "-remove dims from different target area only - zoom" + Opt.figAppend),'file_type',"svg");
poststeps(g)

%% REMOVE PATTERN FROM DIFFERENT TARGET AREA, zoom into [0,6]
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',~rtnew.sameDirection & not_control);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5,'type','sem','type',eb);
% g = g.stat_summary('geom','area','type',eb);
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of different target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval", genH + "-remove dims from different target area only - zoom" + Opt.figAppend),'file_type',"pdf");
g.export('file_name',figuredefine("dimensionRemoval", genH + "-remove dims from different target area only" + Opt.figAppend),       'file_type',"svg");
poststeps(g)

%% Area summary
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'marker', rtnew.targetArea,...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.sameDirectionLabel),'subset',~rtnew.sameDirection & not_control);
g = g.facet_grid([],categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','type',eb);
g = g.stat_summary('geom','area','type',eb);
g = g.set_color_options('chroma',0,'lightness',30);
g.set_text_options("interpreter",'latex');
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Brain area removed');
snapnow;
g=g.update('subset',rtnew.sameDirection & not_control);
%g = g.facet_grid([],categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','type',eb);
g = g.stat_summary('geom','area','type',eb);
g = g.set_color_options(); % Restore default color
set(g.facet_axes_handles, 'xlim', [0 6])
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export("file_name", figuredefine("dimensionRemoval/", genH + "-area-summary" + Opt.figAppend), "file_type", "pdf");
g.export("file_name", figuredefine("dimensionRemoval/", genH + "-area-summary" + Opt.figAppend), "file_type", "svg");
poststeps(g)

    function steps()
        sgtitle(genH);
        set(gcf, 'Position',  get(0, 'Screensize'));
    end
    function poststeps(g)
        % g.draw();
        screensize =  get(0, 'Screensize');
        quarter_screenwidth_size = [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2];
        try
            set(g.fig, 'Position',  quarter_screenwidth_size);
        catch
            set(gcf, 'Position',  quarter_screenwidth_size);
        end
    end

end
