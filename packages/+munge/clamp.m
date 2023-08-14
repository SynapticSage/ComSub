function clamped = clamp(x, xmin, xmax)
    % CLAMP Clamps the values in x between xmin and xmax.
    %
    %   CLAMPED = CLAMP(X, XMIN, XMAX) returns an array of the same size as X 
    %   where values less than XMIN are set to XMIN, values greater than XMAX 
    %   are set to XMAX, and values in between are unchanged.
    %

    if nargin < 3
        error('You must provide x, xmin, and xmax.');
    end

    clamped = max(min(x, xmax), xmin);
end

