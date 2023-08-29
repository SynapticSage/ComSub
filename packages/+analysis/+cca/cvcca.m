function out = cvcca(sp1, sp2, kfold, lambda, maxNpc, sp1_equals_sp2, Npartitions)
    % function out = cvcca(sp1, sp2, kfold, lambda, maxNpc, sp1_equals_sp2, Npartitions)
%
% Inputs: 
% - sp1, sp2: nBins x nNeurons arrays of spike counts
% - kfold: number of folds for cross-validation
% - lambda: regularization parameter for ridge regression
% - maxNpc: maximum number of PCs to use for dimensionality reduction
%
% Outputs: 
% - out: structure with fields:
%   - val: nBins x nBins matrix of correlation values
%   - p: nBins x nBins matrix of p-values
%   - A: nNeurons x maxNpc matrix of canonical vectors for sp1
%   - B: nNeurons x maxNpc matrix of canonical vectors for sp2
%   - cv: cvpartition object used for cross-validation
%   - U: nBins x maxNpc matrix of projections of sp1 onto canonical vectors
%   - V: nBins x maxNpc matrix of projections of sp2 onto canonical vectors
%   - AA: cell array of A matrices for each fold
%   - BB: cell array of B matrices for each fold

if nargin < 7
    Npartitions = 30; % Default value
end

if nargin < 6
    sp1_equals_sp2 = false;
end

% Check for the sp1_equals_sp2 flag
if sp1_equals_sp2

    % Number of neurons
    numNeurons = size(sp1, 2);
    halfN = floor(numNeurons / 2);
    pcs = min(halfN, maxNpc);
    
    AA_accum = cell(Npartitions, 1);
    BB_accum = cell(Npartitions, 1);
    A_accum = zeros(numNeurons, pcs, Npartitions);
    B_accum = zeros(numNeurons, pcs, Npartitions);
    
    for i = progress(1:Npartitions, 'Title', 'Same area partitions')
        % Split into partitions
        idx = randperm(numNeurons);
        new_sp1 = sp1(:, idx(1:halfN));
        new_sp2 = sp1(:, idx(halfN+1:end));

        % Recursive call
        maxNpc1 = min(size(new_sp1, 2), maxNpc);
        out1 = analysis.cca.cvcca(new_sp1, new_sp2, kfold, lambda, maxNpc1, false);
        maxNpc2 = min(size(new_sp2, 2), maxNpc);
        out2 = analysis.cca.cvcca(new_sp2, new_sp1, kfold, lambda, maxNpc2, false);

        if size(out1.A, 2) < pcs 
            pcinds = 1:size(out1.A, 2);
        else
            pcinds = 1:pcs;
        end

        % Placeholder for original dimensions
        A_temp = zeros(numNeurons, pcs);
        B_temp = zeros(numNeurons, pcs);

        % Place the smaller A and B vectors at the respective indices
        A_temp(idx(1:halfN), pcinds)     = out1.A;
        A_temp(idx(halfN+1:end), pcinds) = out2.A;

        B_temp(idx(halfN+1:end), pcinds) = out1.B;
        B_temp(idx(1:halfN), pcinds)     = out2.B;

        A_accum(:,:,i) = A_temp;
        B_accum(:,:,i) = B_temp;

        % % For variable size of AA and BB
        % nCrossVal = numel(out1.AA); % Assuming out1.AA and out2.AA are of the same length
        % for j = 1:nCrossVal
        %     if i == 1
        %         AA_accum{j} = (out1.AA{j} + out2.AA{j})/2;
        %         BB_accum{j} = (out1.BB{j} + out2.BB{j})/2;
        %     else
        %         AA_accum{j} = AA_accum{j} + (out1.AA{j} + out2.AA{j})/2;
        %         BB_accum{j} = BB_accum{j} + (out1.BB{j} + out2.BB{j})/2;
        %     end
        % end
    end

    % % Average over all partitions
    % for j = 1:nCrossVal
    %     AA_accum{j} = AA_accum{j} / Npartitions;
    %     BB_accum{j} = BB_accum{j} / Npartitions;
    % end
    % out.AA = AA_accum;
    % out.BB = BB_accum;
    out.A = mean(A_accum, 3);
    out.B = mean(B_accum, 3);
    out.A = mean(cat(3, out.A, out.B), 3);
    out.B = out.A;

    % Compute U and V for intra-brain interaction
    out.U = sp1 * out.A;
    out.V = sp1 * out.B;

else

    doPCA = true; % if false, just use what you're given

    if nargin<5
        maxNpc = 20;
    end
    if nargin<4
        lambda = [];     
    end
    if nargin<3 || isempty(kfold)
        kfold = 5;
    end

    N = size(sp1,1); %number of trials
    cvp = cvpartition(N, 'KFold', kfold);
    nd = min([maxNpc size(sp2,2), size(sp1,2)]);
    % nd = min([size(sp1,3) size(sp2,3)]);
        
    % reducing dimensionality to help regularize
    X = squeeze(sp1);
    if doPCA        
        [coefX,Xs] = pca(X);
        coefX = coefX(:,1:nd);
        Xs = Xs(:,1:nd);
    else 
        coefX = eye(size(X,2));
        Xs = X;
    end
    Y = squeeze(sp2);
    if doPCA
        [coefY,Ys] = pca(Y);           
        coefY = coefY(:,1:nd);
        Ys = Ys(:,1:nd);
    else
        coefY = eye(size(Y,2));
        Ys = Y;
    end

    % ** New method: use cross-validation
    U = zeros(N,nd);
    V = zeros(N,nd);
    AA = cell(kfold,1);
    BB = cell(kfold,1);

    for k = 1:kfold
        % generate canonical dimension of each matrix on the training set           
        if isempty(lambda)
            [A,B] = canoncorr(...
                Xs(cvp.training(k),:),...
                Ys(cvp.training(k),:));
        else
            % ridge regression version
            [A,B] = canoncorr(...
                [Xs(cvp.training(k),:); lambda*eye(nd); zeros(nd)],...
                [Ys(cvp.training(k),:); zeros(nd); lambda*eye(nd)]);
        end
        if k == 1 && (size(A,2) < nd || size(B,2) < nd)
            nd = min([size(A,2) size(B,2)]);
            U = U(:,1:nd);
            V = V(:,1:nd);
        end
        AA{k} = coefX * A;
        BB{k} = coefY * B;
        
        % project the test set onto the canonical dimension
        U(cvp.test(k),:) = Xs(cvp.test(k),:)*A(:,1:nd);
        V(cvp.test(k),:) = Ys(cvp.test(k),:)*B(:,1:nd);
        
    end

    out.AA = AA;
    out.BB = BB;
    out.cv = cvp;
    out.A = mean(cat(3,AA{:}),3);
    out.B = mean(cat(3,BB{:}),3);
    out.U = U;
    out.V = V;
end

% correlate the projections. since each test set's projections will be zero
% mean already, we can just combine them all here
[rval, pval] = corr(out.U, out.V);
out.val  = rval;
out.pval = pval;

end
