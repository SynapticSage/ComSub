function Patterns = JPECC_RRR(Patterns, Option)
% Computes lagged JPECC between patterns in two brain araes

%Patterns(i).X_source/X_target: [neurons: 63, times: 30, trials]

% Number of cross validation folds.
cvFolds = Option.jpecc.cvnum;

cvOptions = statset('crossval');
regressMethod = @ReducedRankRegress;

nTarget = size(Patterns(1).X_target,1);
nSource = min(size(Patterns(1).X_source,1), size(Patterns(1,2,1).X_source,1));

%numDimsUsedForPrediction = 1:min(nTarget,nSource);
numDimsUsedForPrediction = 
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) ...
    RegressFitAndPredict(regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
    numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
    false, 'Scale', false);


% ------------------------------------------------------------------
% NOTES
% ------------------------------------------------------------------
% RegressMethod is what we use to apply
% RegressFitAndPredict has two lines: regression function and loss eval
% RankRegressRoutine : LIterally applies the crossval function to
% RegressFitAndPredict, and gives the cvLoss matrix and optimal dimension.
% and then recomputes B,V,B_ with optimal dimension
% -------------------------------------------------------------------

% Loop and apply
for i = [3,4,7,8,11,12]
    %progress(1:numel(Patterns), 'Title', 'Rank regression')

    disp(i)
    disp("!!!!!!!!!!!!!!!!!!!!")

    p = Patterns(i);

    tens_curr_source = p.X_source;
    tens_curr_target = p.X_target;


    nBins = size(tens_curr_source,2);

     
    % JPECC : Check all time lags
    for t1 = 1
        disp(t1)
    for t2 = 1
        disp(t2)
        cs = squeeze(tens_curr_source(:, t1, :))';
        ct = squeeze(tens_curr_target(:, t2, :))';

        nan_rows = any(isnan(cs), 2) | any(isnan(ct), 2);

        if(mean(cs, 'all') < 100*eps())
            warning('centering...');
            cs = cs - mean(cs, 2);
            ct = ct - mean(ct, 2);
        end

        cs = cs(~nan_rows, :);
        ct = ct(~nan_rows, :);

        % 计算 s
        s = std(cs,0,1);
        % 计算噪声的幅度
        noise_amplitude = (cvFolds+1)*sqrt(eps(class(s)));
        % 生成与 cs 相同大小的随机噪声
        noise = noise_amplitude * randn(size(cs));
        % 将噪声添加到 cs
        cs = cs + noise;

        lessThan3Samples = sum(~nan_rows) < 3;
        if lessThan3Samples
            p.jpecc(t1, t2).cvl    = nan;
            p.jpecc(t1, t2).cvLoss = nan;
            p.jpecc(t1, t2).optDimReducedRankRegress = nan;
            p.jpecc(t1, t2).B = nan;
            p.jpecc(t1, t2).B_ = nan;
            p.jpecc(t1, t2).V = nan;
            p.jpecc(t1, t2).muOptLoss   = nan;
            p.jpecc(t1, t2).stdOptLoss  = nan;
            p.jpecc(t1, t2).muFullLoss  = nan;
            p.jpecc(t1, t2).stdFullLoss = nan;
            p.jpecc(t1, t2).R2 = nan;
            continue
        end
        
        %{
        % extra condition to avoid error raised when (t1, t2) = (6:30,30) or (30,
        % 6:30)
        if ((t1 >= 6 && t1 <= 30 && t2 == 30) || (t1 == 30 && t2 >= 6 && t2 <= 30))
            p.jpecc(t1, t2).cvl    = nan;
            p.jpecc(t1, t2).cvLoss = nan;
            p.jpecc(t1, t2).optDimReducedRankRegress = nan;
            p.jpecc(t1, t2).B = nan;
            p.jpecc(t1, t2).B_ = nan;
            p.jpecc(t1, t2).V = nan;
            p.jpecc(t1, t2).muOptLoss   = nan;
            p.jpecc(t1, t2).stdOptLoss  = nan;
            p.jpecc(t1, t2).muFullLoss  = nan;
            p.jpecc(t1, t2).stdFullLoss = nan;
            p.jpecc(t1, t2).R2 = nan;
            continue
        end
        %}

        [p.jpecc(t1, t2).cvl,...
         p.jpecc(t1, t2).cvLoss,...
         p.jpecc(t1, t2).optDimReducedRankRegress,...
         p.jpecc(t1, t2).B,...
         p.jpecc(t1, t2).B_,...
         p.jpecc(t1, t2).V] ...
            = rankRegressRoutine(cvFun, cvFolds, ...
            cvOptions, ...
            ct, ...
            cs, ...
            numDimsUsedForPrediction);

         p.jpecc(t1, t2).muOptLoss   = p.jpecc(t1, t2).cvLoss(1, p.jpecc(t1, t2).optDimReducedRankRegress);
         p.jpecc(t1, t2).stdOptLoss  = p.jpecc(t1, t2).cvLoss(2, p.jpecc(t1, t2).optDimReducedRankRegress);
         p.jpecc(t1, t2).muFullLoss  = p.jpecc(t1, t2).cvLoss(1, end);
         p.jpecc(t1, t2).stdFullLoss = p.jpecc(t1, t2).cvLoss(2, end);
         
       
         %calculate explained variance R2
         %{
         B_hat = p.jpecc(t1, t2).B(:, size(ct,2)*(p.jpecc(t1, t2).optDimReducedRankRegress-1)+1:size(ct,2)*p.jpecc(t1, t2).optDimReducedRankRegress);
         [~, Yhat] = RegressPredict(ct, cs, B_hat);

         mean_y = mean(ct, 1);
         SST = sum((ct - mean_y).^2);
         residuals = ct - Yhat;
         SSE = sum(residuals.^2);
         SSR = SST - SSE;

         R2 = SSR / SST;  % 决定系数R-squared
         %}
         [NSE, ~] = RegressPredict(ct,cs,p.jpecc(t1, t2).B);
         p.jpecc(t1, t2).R2 = 1-NSE;

         %{
         B_hat = p.jpecc(t1, t2).B(:, size(ct,2)*(p.jpecc(t1, t2).optDimReducedRankRegress-1)+1:size(ct,2)*p.jpecc(t1, t2).optDimReducedRankRegress);
         [~, Yhat] = RegressPredict(ct, cs, B_hat);
         R = corrcoef(ct, Yhat);
         explained_variance = R(1,2)^2;
         p.jpecc(t1, t2).R2 = explained_variance;
         %}
        
    end
    end

    Patterns(i) = p;

end


