function plot_discrete_colors(colors)
    % colors: Nx3 matrix where each row is an RGB triplet

    num_colors = size(colors, 1);
    
    figure;
    hold on;
    
    for i = 1:num_colors
        rectangle('Position', [i, 0.5, 1, 1], 'FaceColor', colors(i, :), 'EdgeColor', 'none');
    end
    
    xlim([1, num_colors+1]);
    ylim([0, 2]);
    set(gca, 'YTick', []);
    set(gca, 'XTick', 1.5:1:(num_colors+0.5), 'XTickLabel', 1:num_colors);
    title('Discrete Color Palette');
    box on;
end

