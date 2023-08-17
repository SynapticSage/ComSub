
% Calculate correlation coefficients
concatenated_matrix = [selected.pattern.u, selected.pattern.v];
interp_Cavg = interp1(selected.efizz.t, selected.efizz.Cavg, selected.pattern.X_time);
[correlation_matrix, p_values] = corrcoef([interp_Cavg, concatenated_matrix], 'Rows', 'complete');

% Determine significance and set non-significant coefficients to NaN
alpha = 0.05;
non_significant = p_values > alpha;
correlation_matrix(non_significant) = NaN;

% 1. Prepare x-axis labels
num_u = size(selected.pattern.u, 2);
num_v = size(selected.pattern.v, 2);
u_labels = arrayfun(@(x) ['u' num2str(x)], 1:num_u, 'UniformOutput', false);
v_labels = arrayfun(@(x) ['v' num2str(x)], 1:num_v, 'UniformOutput', false);
x_labels = ['Cavg', u_labels, v_labels];

% 2. Prepare y-axis labels
y_labels = arrayfun(@(x) num2str(round(x, 2, 'significant')), efizz.f, 'UniformOutput', false);
y_labels = ['Cavg', y_labels];

% 3. & 4. Plot and customize
imagesc(correlation_matrix);

% Setting labels
set(gca, 'XTick', 1:(num_u + num_v + 1), 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:(length(efizz.f) + 1), 'YTickLabel', y_labels);

% Drawing a white dividing line to separate u and v parts
hold on;
plot([1.5, 1.5], [0.5, length(efizz.f) + 2.5], 'w', 'LineWidth', 2); % Vertical line after Cavg
plot([num_u + 1.5, num_u + 1.5], [0.5, length(efizz.f) + 2.5], 'w', 'LineWidth', 2); % Vertical line separating u and v
plot([0.5, num_u + num_v + 2.5], [1.5, 1.5], 'w', 'LineWidth', 2); % Horizontal line after Cavg
hold off;

% Additional plot customizations
colorbar; % Shows a color scale
title('Correlation Coefficient Matrix');
xlabel('Variables');
ylabel('Frequencies (Hz) & Cavg');

