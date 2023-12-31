function [fig_smoothingthroughout, ...
         fig_smoothingindividual] = frchecks(Spk, varargin)
% function plot_frchecks(Spk)
% Plot the spike rates for each bin to observe time dynamic
%
% INPUTS:
%   Spk: the raw data structure in TheScript.m
%
% OUTPUTS:
%   

disp("Plotting spike rates for each bin to observe time dynamic");
tic

ip = inputParser();
ip.addParameter('saveFigures', false, @islogical);
ip.addParameter('appendFigTitle', '', @ischar);
ip.addParameter('savePath', '', @ischar);
ip.addParameter('visible', 'off', @ischar);
ip.parse(varargin{:});
Opt = ip.Results;

disp("FIG TITLE APPEND = " + Opt.appendFigTitle);

if isempty(Opt.savePath)
    Opt.savePath = fullfile(figuredefine(), "frchecks");
    if ~isfolder(Opt.savePath)
        mkdir(Opt.savePath);
    end
    disp("Save path for figures is " + Opt.savePath);
end

Opt.appendFigTitle = strrep(Opt.appendFigTitle, '_', ' ');
if ~startsWith(Opt.appendFigTitle, ' ')
    Opt.appendFigTitle = [' ' Opt.appendFigTitle];
end

spikeRates = mean(Spk.spikeRateMatrix, 1);
% Assuming your spike rate vector is stored in a variable called 'spikeRates'
windowSize = 10000;
smoothedSpikeRates = smooth(spikeRates, windowSize);
st = Spk.sessionTypePerBin;


% ------------------- Plotting -------------------
% 1st plot shows smoothed rates when smoothing is
% applied throughout the entire session
% -------------------------------------------------
fig_smoothingthroughout = fig("smoothingthroughout" + Opt.appendFigTitle);
set(fig_smoothingthroughout, 'visible', Opt.visible);
tileAxes = tiledlayout(3, 1);
nexttile;
plot(spikeRates, "Color", [1,1,1], "LineWidth", 0.5, "LineStyle", ":")
hold on;
plot(smoothedSpikeRates, "Color", 'b', "LineWidth", 2)
title("Smoothed spike rates" + Opt.appendFigTitle);
nexttile;
plot(smoothedSpikeRates(st == 0));
title("Smoothed spike rates for sleep" + Opt.appendFigTitle);
nexttile;
plot(smoothedSpikeRates(st == 1));
title("Smoothed spike rates for run " + Opt.appendFigTitle);
savefig(fig_smoothingthroughout, fullfile(Opt.savePath, "smoothingthroughout" + Opt.appendFigTitle + ".fig"));
saveas(fig_smoothingthroughout, fullfile(Opt.savePath, "smoothingthroughout" + Opt.appendFigTitle + ".png"));
saveas(fig_smoothingthroughout, fullfile(Opt.savePath, "smoothingthroughout" + Opt.appendFigTitle + ".svg"));


spikeRatesSleep = smoothedSpikeRates(st == 0);
spikeRatesRun = smoothedSpikeRates(st == 1);
spikeRatesSmoothedSleep = smooth(spikeRatesSleep, windowSize);
spikeRatesSmoothedRun = smooth(spikeRatesRun, windowSize);

% ------------------- Plotting -------------------
% 2nd plot shows smoothed rates when smoothing is
% applied individually to sleep and run
% -------------------------------------------------
fig_smoothingindividual =  fig("smoothingindividual" + Opt.appendFigTitle);
set(fig_smoothingindividual, 'Position', get(0, 'Screensize'), 'visible', Opt.visible);
tileAxes = tiledlayout(2, 1);
nexttile;
plot(spikeRatesSleep, "Color", [1,1,1], "LineWidth", 0.5, "LineStyle", ":")
hold on;
plot(spikeRatesSmoothedSleep, "Color", 'b', "LineWidth", 2)
title("Smoothed spike rates for sleep" + Opt.appendFigTitle);
nexttile;
plot(spikeRatesRun, "Color", [1,1,1], "LineWidth", 0.5, "LineStyle", ":")
hold on;
plot(spikeRatesSmoothedRun, "Color", 'b', "LineWidth", 2)
title("Smoothed spike rates for run " + Opt.appendFigTitle);
savefig(fig_smoothingindividual, fullfile(Opt.savePath, "smoothingindividual" + Opt.appendFigTitle + ".fig"));
saveas(fig_smoothingindividual, fullfile(Opt.savePath, "smoothingindividual" + Opt.appendFigTitle + ".png"));
saveas(fig_smoothingindividual, fullfile(Opt.savePath, "smoothingindividual" + Opt.appendFigTitle + ".svg"));

close all;

disp("Plotting spike rates for each bin to observe time dynamic took " + toc + " seconds");

end


