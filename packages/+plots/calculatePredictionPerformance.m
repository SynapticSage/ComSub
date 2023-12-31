function [performance, mu, dev] = calculatePredictionPerformance(X_source, X_target, B)
%
% calculate the predictive performance of source firing to target firing
%
% Input:

[~,nTarget] = size(X_target);
[loss, yhat] = RegressPredict(X_target, X_source, B);
performance = [];
for j = 1:nTarget
    unaccounted_variance = sum((yhat(:,j)-X_target(:,j)).^2);
    total_variance       = sum((mean(X_target(:,j))-X_target(:,j)).^2);
    temp_r_square        = 1-unaccounted_variance/total_variance;

    performance = [performance temp_r_square];
end
mu = mean(performance);
dev = std(performance);

end

