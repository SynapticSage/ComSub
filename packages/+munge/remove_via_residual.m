function A_1_corrected = remove_via_residual(A_1, A_2)
% Remove the influence of A_2 on A_1 by regressing A_1 on A_2 and
% subtracting the predicted values from A_1. This is done for each
% canonical dimension of A_1 separately.


[nNeurons, nCanonicals_A1] = size(A_1);
nCanonicals_A2 = size(A_2, 2);
A_1_corrected = zeros(nNeurons, nCanonicals_A1);

for i = 1:nCanonicals_A1
    y = A_1(:, i);
    % X = [ones(nNeurons, 1), A_2]; % Adding a constant term for the intercept
    X = A_2;
    
    % Estimate regression coefficients
    b = regress(y, X);
    
    % Predicted A_1 values based on A_2
    y_pred = X * b;
    
    % Residuals are the corrected A_1 values after removing influence of A_2
    A_1_corrected(:, i) = y - y_pred;
end

