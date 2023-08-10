%{
-------------------------------
Possible subpsace angle methods
-------------------------------
1. **Principal Angles (Canonical Angles or Subspace Angles):** This method is
probably the most common way to measure the angle between two subspaces.
Principal angles are invariant under unitary transformations. In the context of
linear algebra, the principal angles are defined as the angles between pairs of
nearest points of two subspaces. If we have two subspaces \( U \) and \( V \)
of \( \mathbb{R}^n \) (or \( \mathbb{C}^n \)) with dimensions \( p \) and \( q
\) respectively, there are \( min(p, q) \) principal angles between \( U \) and
\( V \). The cosine of these angles can be found as the singular values of the
cross-covariance matrix \( X^T Y \).
%}


% -------------
% Shortcut vars
% -------------
selectGenH    = "coherence";
directionality = ["hpc-hpc","hpc-pfc"];
patternNames   = ["theta", "delta", "ripples","theta-control",...
                "delta-control", "ripples-control"];
genH           = shortcut.generateH([Option.generateH]);
genH           = unique(genH);
direct         = shortcut.directionality(directionality);
patternSymbols = shortcut.patternSymbols(patternNames, 2);
figureFolder   = fullfile(figuredefine, 'subspaceAngle');

% Reshape so that animals in partition dimension 
% (they function as replicates in this analysis)
[n.animGenH, n.partition, n.dir, n.sing, n.patterns] = size(Patterns);
n.animals = length(unique([Patterns.animal]));
n.genH = length(unique([Patterns.genH_name]));
P=reshape(Patterns, [n.animals,n.genH,n.partition,n.dir,n.patterns]);
assert(length(unique([P(1,:,:,:,:,:).animal])) == 1, "Reshape failed")
disp("Reshape 1: Size of Patterns: " + join(string(size(P)), "x"))
% switch the animal and genH dimensions
P=permute(P,[2,1,3,4,5]);
disp("Permute: Size of Patterns: " + join(string(size(P)), "x"))
% and now finally collapse animals into partition dimension
P=reshape(P,[n.genH,n.animals*n.partition,n.dir,n.patterns]);
disp("Reshape 2: Size of Patterns: " + join(string(size(P)), "x"))

% ----------------------------
% Yank out the relevant data
% ----------------------------
rankRegress = nd.fieldGet(P,'rankRegress');
assert(all(cellfun(@isempty,{rankRegress.B}) == 0), "missing rankRegress.B")
B_ = nd.fieldGetCell(rankRegress,'B_');
V  = nd.fieldGetCell(rankRegress,"V");
sourceIndex = nd.fieldGetCell(P,'index_source');
targetIndex = nd.fieldGetCell(P,'index_target');
smax = max(cellfun(@max, sourceIndex),[],'all');
tmax = max(cellfun(@max, targetIndex),[],'all');


selectGenHs = ["coherence",...
    ..."fft", "wpli"];
    ];
measurements = ["princ",...
        "princ_resamp",...
        "princ_unitvar", ...
        ];

% ----------------------------
% How many components?
% ----------------------------
K = 2:6;
for selectGenH  = progress(selectGenHs,'Title','selectGenH')
for measurement = progress(measurements,'Title','measurement')
    for normalize = [false, true]
        for k = progress(K, 'Title','K')
            % ---------------------
            % Compute dim-reduced B
            % ---------------------
            indices = nd.indicesMatrixForm(rankRegress);
            B = cell(max(indices,[],1));
            for index = progress(indices','Title','Creating B_rrr')
                I = num2cell(index);
                thisB_ = B_{I{:}};
                thisV = V{ I{:} }; % V is the same for all indices
                tmp = thisB_(:,1:k) * thisV(:,1:k)'; % B = B_ * V' (dim-reduced)
                if isempty(selectGenH)
                    tmp = munge.alignToAreaMatrices(tmp, sourceIndex{I{:}}, targetIndex{I{:}}, smax, tmax); %ry 2023
                end
                B{I{:}} = tmp;
            end
            % Munge the shapes used to compute distances here
            partitionDim = 2;
            genHDim = 1; dirDim = 3; patternDim = 4;
            indices = nd.indicesMatrixForm(rankRegress);
            partitions = unique(indices(:,partitionDim));
            % only keep the indices that are in the first 3 partitions
            indices = unique(indices(:,[genHDim,dirDim,patternDim]),'rows');
            indices = indices(indices(:,end)<4,:); % selecting only the first 3 patterns
            assert(~isempty(indices), "No indices found")
            if ~isempty(selectGenH)
                genHs = genH(indices(:, genHDim));
                genHfilter = find(genHs == selectGenH);
                indices = indices(genHfilter, :);
                assert(~isempty(indices), "No indices found")
            end
            subspaceDist = zeros(numel(partitions),...
                            size(indices, 1), size(indices, 1));
            frobDist = zeros(numel(partitions),...
                            size(indices, 1), size(indices, 1));
            rowVar = string([]);
            % expects a row vector ndim == ndims(B) == ndims(P)
            getInd = @(part,rest) [rest(1) part rest(2:end)]; 
            % -----------------------------
            % Compute angle of each subpace
            % -----------------------------
            % (REDO : Once per partition)
            if startsWith(measurement,"princ")
                ang_method = 'principal';
            elseif startsWith(measurement,"frobenius")
                ang_method = 'frobenius';
            else
                error("You spelled it wrong")
            end
                
            assert(~isempty(indices), "No indices found")
            for part = progress(partitions','Title','Partitions')
                for i = 1:size(indices, 1)
                    for j = 1:size(indices, 1)
                        I1 = getInd(part, indices(i,:));
                        I2 = getInd(part, indices(j,:));
                        I1 = num2cell(I1);
                        I2 = num2cell(I2);
                        x = B{I1{:}}; y = B{I2{:}};
                        six = sourceIndex{I1{:}}; tix = targetIndex{I1{:}};
                        siy = sourceIndex{I2{:}}; tiy = targetIndex{I2{:}};
                        assert(isequal(six,siy), "Indices don't match")
                        % subspaceDist(part,i,j) = subspace(x, y, ang_method);
                        [~,r1] = clean.removeDepRows(x);
                        [~,r2] = clean.removeDepRows(y);
                        r = setdiff(1:size(x,1), union(r1,r2));
                        % [~,c1] = clean.removeDepRows(x');
                        % [~,c2] = clean.removeDepRows(y');
                        % c = setdiff(1:size(x,2), union(c1,c2));
                        x = x(r,:); 
                        y = y(r,:);
                        if contains(measurement, "unitvar")
                            % make rows unit vector
                            xr = x ./ vecnorm(x,2,2);
                            yr = y ./ vecnorm(y,2,2);
                            % make columns unit vector
                            xc = x ./ vecnorm(x,2,1);
                            yc = y ./ vecnorm(y,2,1);
                            goodRows = ~any(isnan(xr),2) & ~any(isinf(xr),2) & ~any(isnan(yr),2) & ~any(isinf(yr),2);
                            goodCols = ~any(isnan(xc),1) & ~any(isinf(xc),1) & ~any(isnan(yc),1) & ~any(isinf(yc),1);
                            xr = xr(goodRows,goodCols);
                            xc = xc(goodRows,goodCols);
                            yr = yr(goodRows,goodCols);
                            yc = yc(goodRows,goodCols);
                            subspaceDist(part,i,j) = mean(...
                                [subspaceang(xr, yr, ang_method),...
                                subspaceang(xc, yc, ang_method)]);
                            frobDist(part,i,j) = mean(...
                                [norm(xr-yr,'fro'),...
                                norm(xc-yc,'fro')]);
                        elseif contains(measurement, "resamp")
                            % sample rows and columns 10 times
                            n = 3;
                            tmp = zeros(n,1);
                            tmpf = zeros(n,1);
                            for i = 1:n
                                rows = randsample(1:size(x,1), size(x,1));
                                cols = randsample(1:size(x,2), size(x,2));
                                tmp(i) = subspaceang(x(rows,cols), y(rows,cols), ang_method);
                                tmpf(i) = norm(x(rows,cols)-y(rows,cols),'fro');
                            end
                            subspaceDist(part,i,j) = mean(tmp,1);
                            frobDist(part,i,j)     = mean(tmpf,1);
                        else
                            subspaceDist(part,i,j) = subspaceang(x, y, ang_method);
                            frobDist(part,i,j)     = norm(x-y,'fro');
                        end
                        % figure;tiledlayout(2,1);nexttile; imagesc(x);cmocean('balance');clim([-0.1,0.1]);nexttile;imagesc(y);cmocean('balance');clim([-0.1,0.1]); axis(findobj('type','axes'),'square')
                    end
                    [x, ~, y, z] = deal(I1{:});
                    index     = [direct(y), patternSymbols(z), genH(x)];
                    rowVar(i) = join(string(index),'-');
                    assert(size(subspaceDist,3) == size(subspaceDist,2), "Not square")
                end
            end
            rowVar = cellstr(rowVar);
            disp("Saving subspaceDist with size: " + string(size(subspaceDist)))
            plots.subspace.angle.save(subspaceDist, frobDist, rowVar, k, normalize, measurement, selectGenH);
            tmp = {frobDist,subspaceDist};
            % -----------------------------
            % Checksum of the save file
            % -----------------------------
            [subspaceDist, frobDist, rowVar] = plots.subspace.angle.load(k, normalize, measurement, selectGenH);
            assert(isequal(subspaceDist, tmp{2}), "Subspace dist not equal")
            assert(isequal(frobDist, tmp{1}), "Frobenius dist not equal")
        end
    end
end
end

%% Plot

K = 2:6;
for selectGenH  = selectGenHs
for measurement = measurements
    for normalize = [false, true]
        for k = progress(K, 'Title','K')
            [subspaceDist, frobDist, rowVar] = plots.subspace.angle.load(k, normalize, measurement, selectGenH);
            % -----------------------------
            G = clustergram(subspaceDist, 'ColumnLabels', rowVar, 'RowLabels', rowVar, 'Colormap', crameri('acton'), 'symmetric', true);
%             set(G, 'Colormap', crameri('lajolla'))
%             set(G, 'Colormap', crameri('-grayC'))
            cgFig = findall(0,'Type','Figure','Tag','Clustergram'); %handle to clustergram figure
            cgFig(1).Children(end).YTickLabelRotation=35;
            cgFig(1).Children(end).XTickLabelRotation=-25;
            cgFig(1).Children(end).FontSize=16;
            G.addTitle("subspace distance, K="+k + ", " + measurement + ", " + normalize)
            saveas(cgFig(1),fullfile(figureFolder,sprintf('subspaceDistance_K=%d.svg',k)));
        end % k
    end % normalize
end
end

%ax = findobj(gcf,'type','axes');
%set(ax,'YTickLabelRotation',35, 'XTickLabelRotation', -35, 'FontSize', 20);
%saveas(gcf,fullfile(figureFolder,'subspaceDistance.svg'))
