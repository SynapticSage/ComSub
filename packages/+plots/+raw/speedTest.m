% Extract the required data
theta_time = lfp.hpc.time;
theta_amp = abs(hilbert(lfp.hpc.theta));
behavior_time = behavior.time;
speed = behavior.vel;

% Interpolate theta amplitude to match behavior time points
interp_theta_amp = interp1(theta_time, theta_amp, behavior_time);

% Remove NaNs due to interpolation edges
valid_idx = ~isnan(interp_theta_amp);
speed = speed(valid_idx);
interp_theta_amp = interp_theta_amp(valid_idx);

% Regression
[B_theta,~,~,~,stats_theta] = regress(speed, [ones(size(interp_theta_amp)), interp_theta_amp]);

% Display regression results
disp(['Theta Regression Coefficient: ', num2str(B_theta(2))]);
disp(['Theta R-squared: ', num2str(stats_theta(1))]);

% Extract the required data
efizz_time = efizz.t;
theta_freq_idx = find(efizz.f >= 8, 1); % Assuming 8 Hz is your theta frequency
efizz_theta_amp = abs(hilbert(efizz.S1(:, theta_freq_idx)));

% Interpolate efizz theta amplitude to match behavior time points
interp_efizz_theta_amp = interp1(efizz_time, efizz_theta_amp, behavior_time);

% Remove NaNs due to interpolation edges
valid_idx_efizz = ~isnan(interp_efizz_theta_amp);
speed_efizz = speed(valid_idx_efizz);
interp_efizz_theta_amp = interp_efizz_theta_amp(valid_idx_efizz);

% Regression
[B_efizz,~,~,~,stats_efizz] = regress(speed_efizz, [ones(size(interp_efizz_theta_amp)), interp_efizz_theta_amp]);

% Display regression results
disp(['Theta Efizz Regression Coefficient: ', num2str(B_efizz(2))]);
disp(['Theta Efizz R-squared: ', num2str(stats_efizz(1))]);

% Extract the delta data
delta_time = lfp.hpc.time;
delta_amp = abs(hilbert(lfp.hpc.delta));

% Interpolate delta amplitude to match behavior time points
interp_delta_amp = interp1(delta_time, delta_amp, behavior_time);

% Remove NaNs due to interpolation edges
valid_idx_delta = ~isnan(interp_delta_amp);
speed_delta = speed(valid_idx_delta);
interp_delta_amp = interp_delta_amp(valid_idx_delta);

% Regression
[B_delta,~,~,~,stats_delta] = regress(speed_delta, [ones(size(interp_delta_amp)), interp_delta_amp]);

% Display regression results
disp(['Delta Regression Coefficient for LFP: ', num2str(B_delta(2))]);
disp(['Delta R-squared for LFP: ', num2str(stats_delta(1))]);

% Extract the required indices for the delta range in efizz frequencies
delta_freq_idx = (efizz.f >= 0.5 & efizz.f <= 4);

% Compute the mean of S1 across the specified delta frequency range
efizz_delta_amp = abs(hilbert(mean(efizz.S1(:, delta_freq_idx), 2)));

% Interpolate efizz delta amplitude to match behavior time points
interp_efizz_delta_amp = interp1(efizz_time, efizz_delta_amp, behavior_time);

% Remove NaNs due to interpolation edges
valid_idx_efizz_delta = ~isnan(interp_efizz_delta_amp);
speed_efizz_delta = speed(valid_idx_efizz_delta);
interp_efizz_delta_amp = interp_efizz_delta_amp(valid_idx_efizz_delta);

% Regression
[B_efizz_delta,~,~,~,stats_efizz_delta] = regress(speed_efizz_delta, [ones(size(interp_efizz_delta_amp)), interp_efizz_delta_amp]);

% Display regression results
disp(['Delta Regression Coefficient for Efizz: ', num2str(B_efizz_delta(2))]);
disp(['Delta R-squared for Efizz: ', num2str(stats_efizz_delta(1))]);

% ---- SMOOTHED DELTA ----

% Determine the number of points corresponding to 1 second
sampling_rate_delta = round(1 / mean(diff(lfp.hpc.time)));
window_size_delta = sampling_rate_delta; % 1 second window

% Apply moving average
smooth_delta_amp = movmean(abs(hilbert(lfp.hpc.delta)), window_size_delta);

% Interpolate smoothed amplitude to match behavior time points
interp_smooth_delta_amp = interp1(delta_time, smooth_delta_amp, behavior_time);

% Regression
valid_idx_smooth_delta = ~isnan(interp_smooth_delta_amp);
speed_smooth_delta = speed(valid_idx_smooth_delta);
interp_smooth_delta_amp = interp_smooth_delta_amp(valid_idx_smooth_delta);

[B_smooth_delta,~,~,~,stats_smooth_delta] = regress(speed_smooth_delta, [ones(size(interp_smooth_delta_amp)), interp_smooth_delta_amp]);

% Display regression results
disp(['Smoothed Delta Regression Coefficient for LFP: ', num2str(B_smooth_delta(2))]);
disp(['Smoothed Delta R-squared for LFP: ', num2str(stats_smooth_delta(1))]);

% ---- SMOOTHED THETA ----

% Determine the number of points corresponding to 1 second
sampling_rate_theta = round(1 / mean(diff(lfp.hpc.time)));
window_size_theta = sampling_rate_theta; % 1 second window

% Apply moving average
smooth_theta_amp = movmean(abs(hilbert(lfp.hpc.theta)), window_size_theta);

% Interpolate smoothed amplitude to match behavior time points
interp_smooth_theta_amp = interp1(theta_time, smooth_theta_amp, behavior_time);

% Regression
valid_idx_smooth_theta = ~isnan(interp_smooth_theta_amp);
speed_smooth_theta = speed(valid_idx_smooth_theta);
interp_smooth_theta_amp = interp_smooth_theta_amp(valid_idx_smooth_theta);

[B_smooth_theta,~,~,~,stats_smooth_theta] = regress(speed_smooth_theta, [ones(size(interp_smooth_theta_amp)), interp_smooth_theta_amp]);

% Display regression results
disp(['Smoothed Theta Regression Coefficient for LFP: ', num2str(B_smooth_theta(2))]);
disp(['Smoothed Theta R-squared for LFP: ', num2str(stats_smooth_theta(1))]);


% ---- One from the other ----
% Assuming theta frequency range is already defined as:
theta_freq_range = [6, 12]; % Example: 4-8Hz
[~, theta_min_idx] = min(abs(efizz.f - theta_freq_range(1)));
[~, theta_max_idx] = min(abs(efizz.f - theta_freq_range(2)));

% Calculate the mean theta power in the defined range for efizz.S1
efizz_theta_power = mean(efizz.S1(:, theta_min_idx:theta_max_idx), 2);

% Determine the number of points corresponding to 1 second for lfp.hpc
sampling_rate_theta_hpc = round(1 / mean(diff(lfp.hpc.time)));
window_size_theta_hpc = sampling_rate_theta_hpc; % 1 second window

% Smooth the lfp.hpc theta amplitude using a moving average
smooth_theta_amp_hpc = movmean(abs(hilbert(lfp.hpc.theta)), window_size_theta_hpc);

% Interpolate efizz theta power to match lfp.hpc time points
interp_efizz_theta_power = interp1(efizz.t, efizz_theta_power, lfp.hpc.time);

% Regression: using the interpolated efizz theta as independent variable
% and the smoothed lfp.hpc theta as dependent variable
valid_idx_theta = ~isnan(interp_efizz_theta_power);
smooth_theta_amp_hpc_valid = smooth_theta_amp_hpc(valid_idx_theta);
interp_efizz_theta_power = interp_efizz_theta_power(valid_idx_theta);

[B_theta,~,~,~,stats_theta] = regress(smooth_theta_amp_hpc_valid(:), [ones(size(interp_efizz_theta_power)); interp_efizz_theta_power]');

% Display regression results
disp(['Regression Coefficient between eFizz Theta and LFP Theta: ', num2str(B_theta(2))]);
disp(['R-squared between eFizz Theta and LFP Theta: ', num2str(stats_theta(1))]);

