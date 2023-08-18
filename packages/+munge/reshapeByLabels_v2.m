function Dnew = reshapeByLabels_v2(Dat, dim, labels, varargin)
    ip = inputParser;
    ip.addParameter('checksumSplitField', [], @ischar); % field to checksum, only valid for struct
    ip.parse(varargin{:});
    Opt = ip.Results;

    uniq = struct();
    counts = zeros(1, numel(labels));
    for i = 1:numel(labels)
        uniq.(labels{i}) = unique([Dat.(labels{i})]);
        counts(i) = numel(uniq.(labels{i}));
    end

    % Create an empty version of the struct
    Dnew = nd.emptyLike(Dat(1));
    Dnew = repmat(Dnew, [counts, size(Dat, setdiff(1:ndims(Dat), dim))]);

    % Fill in the new struct
    inds = nd.indicesMatrixForm(Dat); % gets a matrix of indices (rows are complete indices for Dat)
    for ind = inds'
        % Figure out each index for labels
        topinds = zeros(1, numel(labels));
        I = num2cell(ind);
        for i = 1:numel(labels)
            topinds(i) = find(uniq.(labels{i}) == Dat(I{:}).(labels{i}));
        end
        indfull = [topinds(:); ind(setdiff(1:ndims(Dat), dim))];
        indfull = num2cell(indfull);
        Dnew(indfull{:}) = Dat(I{:});
    end
    
end

