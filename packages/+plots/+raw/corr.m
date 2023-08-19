function plot_correlation_matrix(selected, efizz, field, varargin)

ip = inputParser;
ip.addParameter('saveloc', '', @ischar);
ip.addParameter('savetitle', '', @ischar);
ip.parse(varargin{:});
p = ip.Results;

% Calculate correlation coefficients
concatenated_matrix = [selected.pattern.u, selected.pattern.v, selected.pattern.us, selected.pattern.vs];
selected.efizz.(field) = fillmissing(selected.efizz.(field), 'Constant', 0);
interp_Cavg = interp1(selected.efizz.t, selected.efizz.(field), selected.pattern.X_time);
% [correlation_matrix, p_values] = corrcoef([interp_Cavg, concatenated_matrix], 'Rows', 'complete');
[correlation_matrix, p_values] = corr(interp_Cavg, concatenated_matrix, 'Rows', 'complete');

% Determine significance and set non-significant coefficients to NaN
alpha = 0.05;
non_significant = p_values > alpha;
correlation_matrix(non_significant) = NaN;

% 1. Prepare x-axis labels
num_u = size(selected.pattern.u, 2);
num_v = size(selected.pattern.v, 2);
num_us = size(selected.pattern.us, 2);
num_vs = size(selected.pattern.vs, 2);
u_labels = arrayfun(@(x) ['u' num2str(x)], 1:num_u, 'UniformOutput', false);
v_labels = arrayfun(@(x) ['v' num2str(x)], 1:num_v, 'UniformOutput', false);
us_labels = arrayfun(@(x) ['us' num2str(x)], 1:num_us, 'UniformOutput', false);
vs_labels = arrayfun(@(x) ['vs' num2str(x)], 1:num_vs, 'UniformOutput', false);
% x_labels = ['Cavg', u_labels, v_labels];
x_labels = [u_labels, v_labels, us_labels, vs_labels];

% 2. Prepare y-axis labels
y_labels = arrayfun(@(x) num2str(round(x, 2, 'significant')), efizz.f, 'UniformOutput', false);
downs = 10;
y_labels= y_labels(1:downs:end);
% y_labels = ['Cavg', y_labels];

% 3. & 4. Plot and customize
imagesc(correlation_matrix);

% Setting labels
set(gca, 'XTick', 1:(num_u + num_v + num_us + num_vs), 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:downs:(length(efizz.f) + 1), 'YTickLabel', y_labels);

% Drawing a white dividing line to separate u and v parts
hold on;
% plot([0.5, 0.5], [0.5, length(efizz.f) + 2.5], 'w', 'LineWidth', 2); % Vertical line after Cavg
% plot([num_u + 0.5, num_u + 0.5], [0.5, length(efizz.f) + 2.5], 'w', 'LineWidth', 2); % Vertical line separating u and v
% plot([0.5, num_u + num_v + 0.5], [1.5, 1.5], 'w', 'LineWidth', 2); % Horizontal line after Cavg
% Vertical line separating u and v
plot([num_u + 0.5, num_u + 0.5], [0.5, length(efizz.f) + 0.5], 'w', 'LineWidth', 2); 
% Vertical line separating v and us
plot([num_u + num_v + 0.5, num_u + num_v + 0.5], [0.5, length(efizz.f) + 0.5], 'w', 'LineWidth', 2); 
% Vertical line separating us and vs
plot([num_u + num_v + num_us + 0.5, num_u + num_v + num_us + 0.5], [0.5, length(efizz.f) + 0.5], 'w', 'LineWidth', 2); 
hold off;

cmocean('curl');
axis square;

% Additional plot customizations
colorbar; % Shows a color scale
title('Correlation Coefficient Matrix');
xlabel('Variables');
ylabel('Frequencies (Hz)');


if ~isempty(p.saveloc)
    saveas(gcf, fullfile(p.saveloc, p.savetitle + ".fig"))
    saveas(gcf, fullfile(p.saveloc, p.savetitle + ".pdf"))
    saveas(gcf, fullfile(p.saveloc, p.savetitle + ".png"))
    save(fullfile(p.saveloc, p.savetitle + ".mat"), 'correlation_matrix', 'p_values');
end
