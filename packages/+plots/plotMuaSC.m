function plotMuaSC(Spk, varargin)
% plotFR(Spk, varargin)
% plot firing rate of a neuron

% Plots a heatmap of the spike count matrix.
%
% Inputs:
%   Spk           - the results struct
%
% Options:
%   ax          - axis to plot on
%   colorCols   - which color columns to use for the heatmap

ip = inputParser();
ip.addParameter('ax',        gca, @ishandle);
ip.addParameter('colorCols', 1:3, @isnumeric);
ip.addParameter('ylim',      [],  @isnumeric);
ip.addParameter('kws',       {'color','k', 'linestyle', '--'},  @iscell);
ip.addParameter('smooth',    2.5e3,   @isnumeric);
ip.parse(varargin{:});
Opt = ip.Results;

scm   = sum(Spk.spikeCountMatrix, 1);
scm   = scm(Spk.sessionTypePerBin == 1);
times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);

if Opt.smooth > 0
    scm = smooth(scm, Opt.smooth);
end

% Plot the heatmap
axes(Opt.ax);
% hold on;
if isempty(Opt.ylim)
    im=plot(times, scm, Opt.kws{:});
else
    scm = (scm - min(scm))./(max(scm)-min(scm)) .* ...
        (max(Opt.ylim)-min(Opt.ylim));
    im=plot(times, scm, Opt.kws{:});
    ylim(Opt.ylim);
end
uistack(im, 'top');
im.Tag = 'MUA spike count';
% colormap(Opt.ax, cmap);
% caxis([0 1]);
% set(gca, 'YDir', 'normal');
% xlim([times(1) times(end)]);
