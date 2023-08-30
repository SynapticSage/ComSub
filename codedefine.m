function folder = codedefine(varargin)

persistent announced permfolder creation
if ~exist('announced', 'var') || isempty(announced)
    announced = false;
end

% Initialization for creation
if ~exist('creation', 'var')
    creation = false;
end

% Handle the permfolder and creation options
if ~isempty(varargin)
    if varargin{1} == "-permfolder"
        permfolder = varargin(2:end);
        return
    elseif varargin{1} == "-clearpermfolder"
        permfolder = [];
        return
    elseif varargin{1} == "-getpermfolder"
        folder = permfolder;
        return
    elseif varargin{1} == "-creation"
        creation = varargin{2};
        return
    end
end

if ispc
    if ~announced
        disp('Defining code folder for Ziyi machine');
        announced = true;
    end
    folder = 'C:\\Users\BrainMaker\Matlab Drive\Shared\';
elseif ismac
    if ~announced
        disp('Defining code folder for Ryan''s machine');
        announced = true;
    end
    folder = '~/Data/Matlab-DRIVE/Shared/';
elseif isunix
    if ~announced
        disp('Defining code folder for Ryan''s linux machine');
        announced = true;
    end
    folder = '/Volumes/MATLAB-Drive/Shared/';
end

folder = string(folder);
if nargin > 0
    if exist('permfolder', 'var') && ~isempty(permfolder)
        folder = fullfile(folder, permfolder{:}, varargin{:});
    else
        folder = fullfile(folder, varargin{:});
    end
end

% If creation flag is true and directory doesn't exist, create it
if creation
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

end

