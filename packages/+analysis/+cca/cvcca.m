function out = cvcca(sp1, sp2, kfold, lambda, maxNpc)
% function out = cvcca(sp1, sp2, kfold, lambda, maxNpc)
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

% correlate the projections. since each test set's projections will be zero
% mean already, we can just combine them all here
[rval, pval] = corr(U,V);

out.val = rval;
out.pval = pval;
out.A = mean(cat(3,AA{:}),3);
out.B = mean(cat(3,BB{:}),3);
out.cv = cvp;
out.U = U;
out.V = V;
out.AA = AA;
out.BB = BB;

end
