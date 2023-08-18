function Dat = reshapeByLabels(Dat, dim, labels, varargin)
% Supposed we have an a x b x c x d structure, and one of those dimensions dim
% is labeled by labels. This function reshapes the data to be ... x (dim /
% nUniqueLabels) x ...

ip = inputParser;
ip.addParameter('checksumSplitField', [], @ischar); % field to checksum, only valid for struct
ip.parse(varargin{:});
Opt = ip.Results;

% Get the number of unique labels
[uLabels, ia, ic] = unique(labels);

% Assert that the labels have equal counts in ic
assert(all(histcounts(ic, 1:max(ic)+1) == length(ic) / length(uLabels)));
% Assert that all the labels increment by 1
dia = diff(ia);
assert(all(unique(dia) == dia), ...
    'Labels must equally spaced and equal in number');

% Reshape the data
Dat = reshape(Dat, [size(Dat, 1:dim-1), ...
        size(Dat,dim)/length(uLabels),length(uLabels),  ... expanded dims
        size(Dat, dim+1:ndims(Dat))]);


% Checksum the data
if ~isempty(Opt.checksumSplitField)
    disp("Checksumming data " + Opt.checksumSplitField);
    % assert that the checksum field is unique along the dim of interest
    index = repmat({':'}, 1, ndims(Dat));
    index{dim} = 1;
    checksum = [Dat(index{:}).(Opt.checksumSplitField)];
    assert(length(unique(checksum)) == 1,...
        'Checksum field must be unique along the dimension of interest');
    % x=nd.fieldGet(Dat, Opt.checksumSplitField)
    % x(:,:,1,1,1,1,1,1,1)
end

end
