function varargout = load(animal, datatype, varargin)
% General load for all nd branch objs

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('asNd',   false)
ip.addParameter('asTidy', false)
ip.addParameter('get', false);
ip.addParameter('simplefilter', '');
ip.addParameter('level',  []);   % how leveled show the matfile we grab be (0,1,2,3...etc)
ip.addParameter('indices',[]); % used to match indices (day, epoch, tetrode, etc) ... use nan to select all of a level
ip.addParameter('ind',[]); % alias for indices
ip.addParameter('inds',[]); % alias for indices
ip.parse(varargin{:});
Opt = ip.Results;
Opt.animal = animal;

if ~isempty(Opt.ind) && isempty(Opt.indices)
    Opt.indices = Opt.ind;
end
if ~isempty(Opt.inds) && isempty(Opt.indices)
    Opt.indices = Opt.inds;
end

% Whatever  we do in this function, let's make sure we come home
currdir = pwd;
cleanup = onCleanup(@() cd(currdir));

% How  find the files
filefilter = "*" +  Opt.animal + string(datatype) + "*";
folder = ndbFile.folder(animal, datatype, ip.Unmatched);

% Find and sort the files
files = ndbFile.files(string(animal) + string(datatype),...
    Opt.indices, folder,...
    'level', Opt.level);

% Find a super structure
fCount = 0;
if isempty(files)
    warning('No files found');
    res = {};
else
    for file = files'

        % If a simple grep filter is passed
        if ~isempty(Opt.simplefilter) && ~contains(file.name, Opt.simplefilter)
            continue
        end
        if ~endsWith(file.name,'.mat')
            continue
        end

        fCount = fCount + 1;
        if fCount  == 1
            res = load(string([file.folder filesep file.name]));
            res = res.(datatype);
            indices = ndBranch.indicesMatrixForm(res);
            if ~isempty(Opt.indices)
                if ~iscell(Opt.indices)
                    N = min(size(indices,2), size(Opt.indices,2));
                    filt = ismember(indices(:,1:N), Opt.indices(:,1:N), 'rows');
                else
                    N = min(size(indices,2), size(Opt.indices{1},2));
                    filt = cellfun(@(x) ismember(indices(:,1:N), x(:,1:N), 'rows'), Opt.indices, 'UniformOutput', false);
                    filt = cat(2,filt{:});
                    filt = any(filt, 2);
                end
                nullindices = indices(~filt,:);
            else
                nullindices = [];
            end
            for index = nullindices'
                res = ndBranch.set(res, index, []);
            end
        else
            tmp = load(string([file.folder filesep file.name]));
            tmp = tmp.(datatype);
            indices = ndBranch.indicesMatrixForm(tmp);
            if ~isempty(Opt.indices)
                if ~iscell(Opt.indices)
                    N = min(size(indices,2), size(Opt.indices,2));
                    filt = ismember(indices(:,1:N), Opt.indices(:,1:N), 'rows');
                else
                    N = min(size(indices,2), size(Opt.indices{1},2));
                    filt = cellfun(@(x) ismember(indices(:,1:N), x(:,1:N), 'rows'), Opt.indices, 'UniformOutput', false);
                    filt = cat(2,filt{:});
                    filt = any(filt, 2);
                end
                indices = indices(filt,:);
            end
            for index = indices'
                res = ndBranch.set(res, index, ndBranch.get(tmp, index));
            end
        end
    end
end


if Opt.get
    res = ndb.get(res, Opt.indices);
end
if Opt.asNd && Opt.asTidy
    error('Not implemented yet. Will be ndb -> nd -> tidy');
elseif Opt.asNd 
    res = ndBranch.toNd(res);
elseif Opt.asTidy
    Params = ip.Unmatched;
    Params.datatype = datatype;
    res = tidyData.fromNdBranch_tidy(res, Params);
end

% Output or load caller workspace
if nargout == 0
    assignin('caller',datatype,res);
else
    varargout{1} = res;
    if nargout == 2 && Opt.asNd == false
        varargout{2} = ndBranch.type(res);
    end
end
