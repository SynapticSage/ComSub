function fill_curve(X, Y, varargin)
    % If not specified, default color is blue

    if iscolumn(X)
        X = X';
    end
    if iscolumn(Y) && isvector(Y)
        Y = Y';
    end
    if ismatrix(Y) && size(Y,1) > size(Y,2)
        Y = Y'; 
        [Y, lower] = deal(Y(1,:), Y(2,:));
    else
        Y = Y(1,:);
        lower = zeros(1, length(Y));
    end

    % if nargin < 3
    %     fill_color = 'b';
    % end

    % Create patch vertices
    x_patch = [X, fliplr(X)]; % X values for patch (tracing the curve forward then backward)
    y_patch = [Y, lower]; % Y values for patch (curve values, then back down to y=0)

    % Create the filled patch
    patch(x_patch, y_patch, varargin{:});

    % % Overlay the curve on top (optional)
    % hold on;
    % plot(X, Y, 'k-');
    % hold off;
end

