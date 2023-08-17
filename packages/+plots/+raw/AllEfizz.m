% SEE toplevel script plots.event.plot2016

f = fig(animal + " 2016 Figure epoch " + epoch_option + " " + sortprop + " trajbound " + trajbound_option);
figdict = containers.Map();
figdict("figure") = f;
clf;

% Define the relative heights
rel_heights = [2, 1, 0.5, 0.5, 0.5, 0.5];
names = ["HPC", "PFC", "Theta", "Ripple","U", "V"];
normalized_heights = rel_heights / sum(rel_heights);
cumulative_heights = [0, cumsum(normalized_heights)];
% Define user-defined offsets (if needed)
user_offset.fft.theta.hpc   = -5;    
user_offset.fft.theta.pfc   = 5;   
user_offset.fft.theta.wpli  = 0;
user_offset.fft.theta.cavg  = 0;
user_offset.fft.ripple.hpc  = -7.5; 
user_offset.fft.ripple.pfc  = 7.5;  
user_offset.fft.ripple.wpli = 0;    
user_offset.fft.ripple.cavg = 0;    
user_offset.lfp.theta.hpc   = -5;   
user_offset.lfp.theta.pfc   = 5;    
user_offset.lfp.theta.wpli  = 0;    
user_offset.lfp.theta.cavg  = 0;    
user_offset.lfp.ripple.hpc  = -7.5; 
user_offset.lfp.ripple.pfc  = 7.5;  
user_offset.lfp.ripple.wpli = 0;    
user_offset.lfp.ripple.cavg = 0;    
% Compute auto offsets based on the middle of the visualization range
visualize.theta    = [0, 40];
visualize.ripple   = [130, max(efizz.f)];
auto_offset.theta  = mean(visualize.theta);
auto_offset.ripple = mean(visualize.ripple);

ax = gobjects(length(rel_heights), 1); % Pre-allocate memory for axis handles
rowcnt = 0;
for i = 1:length(rel_heights)
    % Create a subplot based on calculated positions
    if rel_heights(i) > 0
        rowcnt = rowcnt + 1;
        figdict(names(i)) = rowcnt;
    else
        continue;
    end
    disp("plotting " + names(i) + " on row " + rowcnt + ", plottype = " + i);
    figure(figdict("figure"));
    ax(rowcnt) = subplot('Position', [0.1, cumulative_heights(rowcnt) + 0.01, 0.8, normalized_heights(i) - 0.02]);

    switch char(names(i))
        case 'HPC'
            % Spike Raster for Hippocampal Cells
            % figure;tiledlayout('flow');nexttile;imagesc(selected.pattern.a);nexttile;imagesc(sel);
            cla
            spk = selected.spike.times_spiking(cells.hpc);
            [~, isort] = sort(sortby(cells.hpc), sortdir);
            if colorbycomp
                if pattern_overall_ind ~= numel(Patterns_overall)
                    error("colorbycomp only works for pattern_overall_ind = " + numel(Patterns_overall));
                end
                a = selected.pattern.a;
                stds = std(a, [], 1);
                sel = bsxfun(@gt, a, 1*stds) | bsxfun(@lt, a, -1*stds);
                assert(size(sel,1) == size(spk,2));
                for j = 1:size(spk,2)
                    if isempty(spk{j})
                        continue;
                    end
                    sptmp = spk;
                    [sptmp{1:j-1}] = deal([nan]);
                    [sptmp{j+1:end}] = deal([nan]);
                    color = [0 0; 0 1; 1 0] * sel(j,1:2)';
                    darken = 1.2;
                    color = color ./ (norm(color)*darken);
                    color = fillmissing(color, 'constant', 0);
                    lineformat = struct('Color', color);
                    hold on
                    plotSpikeRaster(sptmp, 'PlotType', 'vertline', 'LineFormat', lineformat, ...
                        'VertSpikePosition', isort(j));
                    hold on;
                end
            else
                plotSpikeRaster(spk(isort), 'PlotType', 'vertline', 'AutoLabel', false);
                get(gca,'Children');
            end
            ylabel('HPC Neuron');
            title('Spike Raster of Hippocampal Cells');
        case 'PFC'
            % Spike Raster for Prefrontal Cells
            spk = selected.spike.times_spiking(cells.pfc);
            [~, isort] = sort(sortby(cells.pfc), sortdir);
            if colorbycomp
                if pattern_overall_ind ~= numel(Patterns_overall)
                    error("colorbycomp only works for pattern_overall_ind = " + numel(Patterns_overall));
                end
                b = selected.pattern.b; % FIXED on 2023-08
                stds = std(b, [], 1);
                sel = bsxfun(@gt, b, 1*stds) | bsxfun(@lt, b, -1*stds);
                assert(size(spk,2) == size(sel,1));
                for j = 1:size(spk,2)
                    if isempty(spk{j})
                        continue;
                    end
                    sptmp = spk;
                    [sptmp{1:j-1}] = deal(nan);
                    [sptmp{j+1:end}] = deal(nan);
                    color = [0 0; 0 1; 1 0] * sel(j,1:2)';
                    darken = 1.2;
                    color = color ./ (norm(color)*darken);
                    color = fillmissing(color, 'constant', 0);
                    lineformat = struct('Color', color);
                    hold on
                    plotSpikeRaster(sptmp, 'PlotType', 'vertline', 'LineFormat', lineformat, ...
                        'VertSpikePosition', isort(j));
                    hold on;
                end
            else
                plotSpikeRaster(spk(isort), 'PlotType', 'vertline', 'AutoLabel', false);
                get(gca,'Children');
            end
            ylabel('PFC Neuron');
            title('Spike Raster of Prefrontal Cells');
        case 'Theta'
            tmp = log10(abs(selected.efizz.S1'));
            imagesc(selected.efizz.t, efizz.f, tmp); 
            set(gca, 'YDir', 'normal', 'YLim', visualize.theta); 
            [M,m] = deal(quantile(tmp(:, theta_idx),0.97,'all'), quantile(tmp(:, theta_idx),0.05,'all'));
            cmocean('thermal');
            clim([m, M]);
            hold on;
            if use_fft_avg
                scale = median(structfun(@abs,user_offset.theta));
                set(gca, 'Layer', 'top');
                plot(selected.efizz.t, scale*selected.efizz.theta_hpc + auto_offset.theta + user_offset.fft.theta.hpc, 'w', 'DisplayName', 'HPC Theta', 'LineWidth', 1.5);
                set(gca, 'Layer', 'top');
                plot(selected.efizz.t, scale*selected.efizz.theta_pfc + auto_offset.theta + user_offset.fft.theta.pfc, 'r', 'DisplayName', 'PFC Theta', 'LineWidth', 1.5);
                plot(selected.efizz.t, scale*selected.efizz.theta_wpli + auto_offset.theta + user_offset.fft.theta.wpli, 'k', 'DisplayName', 'WPLI Theta', 'LineWidth', 1.5, 'LineStyle', ':');
                plot(selected.efizz.t, scale*selected.efizz.theta_cavg + auto_offset.theta + user_offset.fft.theta.cavg, 'k', 'DisplayName', 'C-Avg Theta', 'LineWidth', 1.5, 'LineStyle', '--');
            else
                scale = median(structfun(@abs,user_offset.theta));
                plot(selected.lfp.time, scale*selected.lfp.hpc.theta + auto_offset.theta + user_offset.lfp.theta.hpc, 'DisplayName', 'HPC Theta', 'LineWidth', 1.5, 'Color', [1 1 1 0.5]);
                plot(selected.lfp.time, scale*selected.lfp.pfc.theta + auto_offset.theta + user_offset.lfp.theta.pfc, 'DisplayName', 'PFC Theta', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]);
                plot(selected.efizz.t, scale*selected.efizz.theta_wpli + auto_offset.theta + user_offset.fft.theta.wpli, 'k', 'DisplayName', 'WPLI Theta', 'LineWidth', 1.5, 'LineStyle', ':');
                plot(selected.efizz.t, scale*selected.efizz.theta_cavg + auto_offset.theta + user_offset.fft.theta.cavg, 'k', 'DisplayName', 'C-Avg Theta', 'LineWidth', 1.5, 'LineStyle', '--');
            end
            ylabel('Frequency (Hz)');
            title('Theta Band Average Power for HPC and PFC');
            legend('Location', 'northwest');
        case 'Ripple'
            tmp = log10(abs(selected.efizz.S1'));
            imagesc(selected.efizz.t, efizz.f, tmp);
            [M,m] = deal(quantile(tmp(:, ripple_idx),0.7,'all'), quantile(tmp(:, ripple_idx),0.01,'all'));
            cmocean('thermal');
            clim([m, M]);
            set(gca, 'YDir', 'normal', 'YLim', visualize.ripple);
            hold on;
            if use_fft_avg
                scale = 2*median(structfun(@abs,user_offset.ripple));
                plot(selected.efizz.t, scale*selected.efizz.ripple_hpc + auto_offset.ripple + user_offset.fft.ripple.hpc, 'DisplayName', 'HPC Ripple', 'LineWidth', 1.5, 'Color', [1 1 1 0.5]);
                plot(selected.efizz.t, scale*selected.efizz.ripple_pfc + auto_offset.ripple + user_offset.fft.ripple.pfc, 'DisplayName', 'PFC Ripple', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]);
            else
                scale = median(structfun(@abs,user_offset.ripple));
                plot(selected.lfp.time, scale*selected.lfp.hpc.ripple + auto_offset.ripple + user_offset.lfp.ripple.hpc, 'w', 'DisplayName', 'HPC Ripple', 'LineWidth', 1.5, 'Color', [1 1 1 0.5]);
                plot(selected.lfp.time, scale*selected.lfp.pfc.ripple + auto_offset.ripple + user_offset.lfp.ripple.pfc, 'r', 'DisplayName', 'PFC Ripple', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]);
                plot(selected.efizz.t, scale*selected.efizz.ripple_wpli + auto_offset.ripple + user_offset.fft.ripple.wpli, 'k', 'DisplayName', 'WPLI Ripple', 'LineWidth', 1.5, 'LineStyle', ':');
                plot(selected.efizz.t, scale*selected.efizz.ripple_cavg + auto_offset.ripple + user_offset.fft.ripple.cavg, 'k', 'DisplayName', 'C-Avg Ripple', 'LineWidth', 1.5, 'LineStyle', '--');
            end
            ylabel('Frequency (Hz)');
            title('Ripple Band Average Power for HPC and PFC');
        case 'U'
            % Shading for U Component 1
            % fill(X, Y1, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            hold on;
            plot(selected.pattern.X_time, selected.pattern.u(:,1), 'DisplayName', 'HPC U');
            plots.fill_curve(selected.pattern.X_time, selected.pattern.u(:,1), 'b')
            clear alpha
            alpha(0.5);

            % Shading for U Component 2
            % fill(X, Y2, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(selected.pattern.X_time, selected.pattern.u(:,2), 'DisplayName', 'HPC U, 2', 'Color', [0 1 0 0.35], 'LineWidth', 0.5);
            plots.fill_curve(selected.pattern.X_time, selected.pattern.u(:,2), 'g')
            alpha(0.3);
            if dim3
                % Shading for U Component 3
                % fill(X, Y3, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                plot(selected.pattern.X_time, selected.pattern.u(:,dim3), 'DisplayName', 'HPC U, 3', 'Color', [1 0 0 0.2], 'LineWidth', 0.5);
                plots.fill_curve(selected.pattern.X_time, selected.pattern.u(:,3), 'r')
                alpha(0.1);
            end
            sumabs = true;
            if sumabs
                % add in light dotted black the sum of hte absolute values of the top 3 components
                plot(selected.pattern.X_time, sum(abs(selected.pattern.u(:,1:3)),2), 'DisplayName', 'HPC U, 1-3', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'LineStyle', ':');
            end
            yline(0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
            ylabel('U Component');
            title('U Component of hpc-pfc communication');

        case 'V'
            % Shading for V Component 1
            plot(selected.pattern.X_time, selected.pattern.v(:,1), 'DisplayName', 'PFC V');
            hold on;
            plots.fill_curve(selected.pattern.X_time, selected.pattern.v(:,1), 'b');
            alpha(0.5);
            if dim3
                % Shading for V Component 3
                plot(selected.pattern.X_time, selected.pattern.v(:,dim3), 'DisplayName', 'PFC V, 3', 'Color', [1 0 0 0.2], 'LineWidth', 0.5);
                plots.fill_curve(selected.pattern.X_time, selected.pattern.v(:,3), 'r');
                alpha(0.1);
            end
            sumabs = true;
            if sumabs
                % add in light dotted black the sum of hte absolute values of the top 3 components
                plot(selected.pattern.X_time, sum(abs(selected.pattern.v(:,1:3)),2), 'DisplayName', 'PFC V, 1-3', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'LineStyle', ':');
            end
            % Shading for V Component 2
            yline(0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
            plot(selected.pattern.X_time, selected.pattern.v(:,2), 'DisplayName', 'PFC V, 2', 'Color', [0 1 0 0.35], 'LineWidth', 0.5);
            plots.fill_curve(selected.pattern.X_time, selected.pattern.v(:,2), 'g');
            clear alpha
            alpha(0.3);
            ylabel('V Component');
            title('V Component of hpc-pfc communication');

        case 'Raw'

            plot(lfp.hpc.time, lfp.hpc.data(:,1));
            hold on
            plot(lfp.pfc.time, lfp.pfc.data(:,1));
            title('Raw Voltage');
            ylabel('Amplitude');
        
        case 'RawTheta'

            plot(lfp.hpc.theta.time, lfp.hpc.theta.data(:,1));
            hold on
            plot(lfp.pfc.theta.time, lfp.pfc.theta.data(:,1));
            title('Theta Power');
            ylabel('Amplitude');

        case 'RawRipple'

            plot(lfp.hpc.theta.time, lfp.hpc.theta.data(:,1));
            plot(lfp.pfc.theta.time, lfp.pfc.theta.data(:,1));
            title('Ripple Power');
            ylabel('Amplitude');

    end

    % Network pattern
    if ~isempty(shadeOption)
        % Get the current axis
        thisax = gca;
       % Make the current axis active
        hold on;
        % Iterate over the windows in cellOfWindows
        for netpat = 1:numel(selected.cellOfWindows)
        if ismember(netpat, wins)
        for win_idx = 1:size(selected.cellOfWindows{netpat}, 1)
            currWin = selected.cellOfWindows{netpat}(win_idx, :);
            % Shade the region using patch
            if colorbycomp
                color = compcolors(netpat, :);
            else
                color = [0.8, 0.8, 0.8];
            end
            yLimits = ylim(thisax);
            hPatch = patch('XData', [currWin(1), currWin(1), currWin(2), currWin(2)], ...
                  'YData', [yLimits(1), yLimits(2), yLimits(2), yLimits(1)], ...
                  'FaceColor', color, 'EdgeColor', 'none', ...
                  'FaceAlpha', 0.3, 'Parent', ax(rowcnt));
                        % Ensure this patch is not added to legend
            hAnnotation = get(hPatch, 'Annotation');
            legendInfo = get(hAnnotation, 'LegendInformation');
            set(legendInfo, 'IconDisplayStyle', 'off');
            hold on;
        end
        end
        end
        % Restore the hold state of the axis
        hold off; 
    end

end
linkaxes(ax, 'x');  % Link all the x-axes together

% Overlay Linear Distance on the HPC spike raster
lindist_to_sections = ["HPC", "U", "V"];
for sect = lindist_to_sections
    if figdict.isKey(sect) % Check if HPC subplot exists
        % figure(figdict(sect));
        ax2 = axes('Position', get(ax(figdict(sect)), 'Position'), 'Color', 'none');
        yyaxis(ax2, 'right');
        hold on;
        % Define the limits for the secondary Y-axis (adjust these as per your needs)
        y2_lim = [min(behavior.lindist), max(behavior.lindist)];
        ylim(ax2, y2_lim);
        p1 = plot(selected.behavior_time, selected_rows.lindist,...
        'Color', [0 0 0 0.35], 'Parent', ax2, 'LineWidth', 3.5);
        % set(p1, 'Color', [1 0 0 0.5]);  % Here, [1 0 0] is red color and 0.35 is the alpha (transparency)
        ylabel('Linear Distance');
        % ax2.YAxisLocation = 'right';  % Ensure the secondary y-axis is on the right
        ax2.Box = 'off';  % This removes the box to make it cleaner
        linkaxes([ax(figdict("HPC")), ax2], 'x');  % Link the x-axes together
    end
end
sgt = gcf;
set(gcf, 'Position', get(0, 'Screensize'));  % Maximize the figure window
% plots.positionFigOnRightMonitor(gcf)
r = @(x) string(replace(replace(replace(x, " ", "_"), ":", "_"), newline, "_"));
saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.png'));
saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.pdf'));
saveas(gcf, fullfile(savefolder, r(sgt.Name) + '.fig'));

% width = 600;
% height = 1200;
% position = [100, 100, width, height];
% set(gcf, 'Position', position);
