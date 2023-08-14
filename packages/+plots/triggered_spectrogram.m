function triggered_spectrogram(efizz, specs, oneDplots, varargin)
% plots triggered spectrogram data from efizz

ip = inputParser;
ip.addParameter('time',       [], @isnumeric);
ip.addParameter('thresholds', [], @isstruct);
ip.addParameter('means',      [], @isstruct);
ip.addParameter('zscore',     [], @(x) iscellstr(x) || isstring(x));
ip.addParameter('nolog',      [], @(x) iscellstr(x) || isstring(x));
ip.addParameter('freqs',      [], @isnumeric); % frequencies to mean over if plotting specs
ip.addParameter('subtract_mins', true, @islogical);
ip.addParameter('log_yaxis', false, @islogical);
ip.addParameter('freq_ylims',      [], @isnumeric);
ip.addParameter('upperylim', 160, @isnumeric);
ip.parse(varargin{:});
Opt = ip.Results;
linestyles = ["-", "--", "-.", ":"];

% figure
fields = fieldnames(specs);
if isempty(Opt.time)
    Opt.time = 1:size(specs.(fields{1}),1);
end
tiledlayout(length(fields)+length(oneDplots), 1)
fg = gcf; clf;
% set(fg, 'Position', get(0,'Screensize'))
axs = gobjects(length(fields), 1);
for i = 1:length(fields)
    axs(i) = nexttile;
    f = fields{i};
    s = specs.(f);
    if ~isempty(Opt.means) && isfield(Opt.means, f) && f ~= "phi"
        s = s - Opt.means.(f);
    end
    if Opt.subtract_mins && f ~= "phi"
        s = s - min(s, [], 1);
    end
    Opt.time = Opt.time - mean(Opt.time); % center time axis
    if f == "wpli_avg" || f == "Cavg"
        s = imgaussfilt(s, 2);
    end
    if ismember(f, Opt.zscore)
        % along the time axis
        specs.(f) = zscore(s, [], 1);
        imagesc(Opt.time, efizz.f, s')
    elseif ismember(f, Opt.nolog)
        imagesc(Opt.time, efizz.f, s')
    else
        imagesc(Opt.time, efizz.f, signedlog(s)')
    end
    if f == "phi"
        clim([-pi pi])
        cmocean('phase')
        colorbar;
    end
    set(gca, 'YDir', 'normal')
    if ~isempty(Opt.upperylim)
        set(gca, 'YLim', [min(ylim()) Opt.upperylim])
    end
    title(f)
    if ~isempty(Opt.freqs)
        fcnt = 0;
        for freq = Opt.freqs
            fcnt = fcnt + 1;
            y = mean(s(:, efizz.f >= freq(1) & efizz.f <= freq(2)), 2);
            hold on;
            plot(Opt.time, y, 'Color', 'white', 'LineWidth', 2, 'LineStyle', linestyles(mod(fcnt-1, length(linestyles))+1))
        end
    end
end
linkaxes(axs, 'xy')
if ~isempty(Opt.freq_ylims)
    set(gca, 'YLim', Opt.freq_ylims)
end
if Opt.log_yaxis
    set(findobj(gcf,'type','axes'), 'YScale', 'log')
end
if ~isempty(oneDplots)
    for i = 1:length(oneDplots)
        ax=nexttile;
        if iscell(oneDplots)
            odp = oneDplots{i};
        else
            odp = oneDplots(i);
        end
        fields = fieldnames(odp);
        for j = 1:length(fields)
            if ismember(fields{j}, Opt.zscore)
                % along the time axis
                odp.(fields{j}) = zscore(odp.(fields{j}), [], 1);
            end
            if ~isempty(Opt.thresholds) && isfield(Opt.thresholds, fields{j})
                yline(Opt.thresholds.(fields{j}), 'r')
            end
            time = linspace(min(Opt.time), max(Opt.time), length(odp.(fields{j})));
            plot(time, odp.(fields{j}))
            hold on
        end
        set(gca,'xlim', [min(Opt.time) max(Opt.time)])
        legend(fields)
        title(strjoin(fields, ' '))
        linkaxes([axs(end-1) ax], 'x')
    end
end
