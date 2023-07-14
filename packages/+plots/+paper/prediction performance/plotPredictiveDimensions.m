function [full_model, pred_by_dim] = plotPredictiveDimensions(numDimsUsedForPrediction, cvLoss, varargin)

% Plot Reduced Rank Regression cross-validation results or factor analysis
% cross-validation results
% also returns the averaged full model
% full_model = 1-curr_cvLoss(1,end);


ip = inputParser;
ip.addParameter('mode', "rr", @isstring);
ip.addParameter('color', "black", @isstring);
ip.addParameter('averaged', true, @islogical);
ip.addParameter('normalized', false, @islogical);
ip.addParameter("do_plot", true, @islogical);
ip.addParameter("optDim", []);
ip.parse(varargin{:});
opt = ip.Results;



x = 1:numDimsUsedForPrediction;
nonsingular = [];
if ~opt.averaged
    if iscell(cvLoss)
        cvLoss = cell2mat(cvLoss);
    end
    y = 1-cvLoss(1,:);
    e = cvLoss(2,:);
    full_model = 1-cvLoss(1,end);
    if ~isempty(opt.optDim)
        optDim = opt.optDim;
    end
   
    
    for i = 1:size(cvLoss,2)
        pred_by_dim(i) = 1-cvLoss(1,i);
    end
else % average the performance for that partition
    numPartitions = numel(cvLoss);
    toAverage = [];
    
    optDim = [];
    error = [];
    all_full_model = [];
    
    for i = 1:numPartitions
        if ~isempty(opt.optDim)
            curr_optDim = opt.optDim(i);
            optDim = [optDim, curr_optDim];
        end
        curr_cvLoss = cvLoss{i};
        if ~isempty(curr_cvLoss)
            nonsingular = [nonsingular,i];
        end
    end
    optDim = min(cell2mat(optDim));
    
    for i = 1:numel(nonsingular)
        curr_cvLoss = cvLoss{nonsingular(i)};
        try
            toAverage = [toAverage; curr_cvLoss(1,1:optDim)];
        catch
            keyboard
        end
        pred_by_dim (i,:) = 1-curr_cvLoss(1,:);
        error = [error; curr_cvLoss(2,1:optDim)];
        all_full_model = [all_full_model, 1-curr_cvLoss(1,end)];
    end
    y = 1-mean(toAverage,1);
    e = mean(error,1);
    pred_by_dim = mean(pred_by_dim,1);
    full_model = mean(all_full_model);
end

x = 1:numel(y);
if opt.mode == "fa"
    x(x > optDim) = [];
end

if opt.normalized
    y_max = y(end);
    
    y = y./y_max;
    e = e./y_max;
end

%% if asking for a plot, produce it
if opt.do_plot
    errorbarplot = errorbar(x, y, e);
    %     errorbarT(errorbarplot, .5, 2);
    errorbarplot.Color = opt.color;
    
    xlabel('# predictive dimensions')
    
    if opt.normalized
        ylabel("norm. performance")
    else
        ylabel('Performance')
    end
    
    % plot where the optimal dimension falls
    if ~isempty(opt.optDim)
        if ~opt.normalized
            lineObject=line([optDim,optDim],[0 full_model]);
        else
            lineObject =line([optDim,optDim],[0 1]);
        end
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 2;  % Thicken the line
        lineObject.Color = opt.color; % Color it
    end
    
    hold on
    
    plot(1, full_model,'^');
    
end

end