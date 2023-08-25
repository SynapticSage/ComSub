load(Option.animal + "spectralBehavior.mat");
if isfield(efizz, "Cavg")
    disp("Using new average coherence field");
    spectrogram   = efizz.Cavg;
else
    disp("Using old coherence field");
    spectrogram   = efizz.C;
end
frequencyAxis = efizz.f;
times         = efizz.t;
[H, Hvals, Hnanlocs, times] = events.generateFromSpectra(times, spectrogram, frequencyAxis,...
    Option.frequenciesPerPattern);
if Option.sourceArea == "CA1"
    spectrogram = efizz.S1;
else
    spectrogram = efizz.S2;
end
Cavg = efizz.Cavg;
[runningSessions, sleepSessions] = getRunningSessions(Option.animal);
[Hpower, ~, ~, ~] = events.generateFromSpectra(times, spectrogram, ...
    frequencyAxis, Option.frequenciesPerPattern)
    
[Hcoherence, Hvals, Hnanlocs, times] = events.generateFromSpectra(times, Cavg, frequencyAxis,...
    Option.frequenciesPerPattern);
% Get the median power
medianPower = median(Hpower, 'all', 'omitnan');
% Get the coherence
H = Hcoherence;
% Replace with NaN where the power is less than the median power
H(Hpower > medianPower)        = NaN;
Hvals(Hpower > medianPower)    = 0;
Hnanlocs(Hpower > medianPower) = NaN;
theta_indices = find(frequencyAxis >= 6 & frequencyAxis <= 12);
theta_mean_S1 = mean(spectrogram(:, theta_indices), 2);
theta_mean_Cavg = mean(Cavg(:, theta_indices), 2);

% 2D histogram
figure;
histogram2(theta_mean_S1, theta_mean_Cavg, 'DisplayStyle','tile','ShowEmptyBins','on');
colorbar;
title('2D Histogram of Mean Theta Power in S1 vs Cavg');
xlabel('Mean Theta Power in S1');
ylabel('Mean Theta Power in Cavg');

% Calculate the quantiles for each variable
quantiles_S1 = quantile(theta_mean_S1, [0.15, 0.5, 0.85]);
quantiles_Cavg = quantile(theta_mean_Cavg, [0.15, 0.5, 0.85]);

% Add white vertical lines at the quantile values for theta_mean_S1
for i = 1:length(quantiles_S1)
    line([quantiles_S1(i), quantiles_S1(i)], ylim, 'Color', 'white', 'LineStyle', '--');
end

% Add white horizontal lines at the quantile values for theta_mean_Cavg
for i = 1:length(quantiles_Cavg)
    line(xlim, [quantiles_Cavg(i), quantiles_Cavg(i)], 'Color', 'white', 'LineStyle', '--');
end

% Check for directory existence and create if necessary
if ~exist(figuredefine("spec_2dhist", ""), 'dir')
    mkdir(figuredefine("spec_2dhist", ""));
end

saveas(gcf, figuredefine("spec_2dhist", Option.animal + "thetaPowerS1vsCavg.png"))
saveas(gcf, figuredefine("spec_2dhist", Option.animal + "thetaPowerS1vsCavg.fig"))


[~, rank_S1] = sort(theta_mean_S1);
[~, rank_Cavg] = sort(theta_mean_Cavg);
num_samples = length(theta_mean_S1);  % assuming S1 and Cavg are of the same length
norm_rank_S1 = rank_S1 / num_samples;
norm_rank_Cavg = rank_Cavg / num_samples;


figure;
histogram2(norm_rank_S1, norm_rank_Cavg, 'DisplayStyle','tile','ShowEmptyBins','on', 'NumBins', 1000);
colorbar;
title('2D Histogram of Normalized Ranks for Theta Power in S1 vs Cavg');
xlabel('Normalized Rank for Theta Power in S1');
ylabel('Normalized Rank for Theta Power in Cavg');



