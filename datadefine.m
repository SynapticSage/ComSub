function folder = datadefine(varargin)
% Automates pointing to the data  folders location on our computers

disp("Running datadefine")

if ispc
    disp('Defining data folder for Ziyi machine');
    folder = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
elseif ismac
    disp('Defining data folder for Ryan''s machine');
    %folder = '~/Data/commsubspace/';
    folder(1) = "/Volumes/MATLAB-Drive/Shared/SingleDayExpt/"
    folder(2) = "~/Data/commsubspace";
elseif isunix
    disp('Defining data folder for Ryan''s linux machine');
    folder = fullfile(codedefine, 'SingleDayExpt');
end

if ~isempty(varargin)
    folder = fullfile(folder, varargin{:});
end

disp("Folders =")
disp(folder)

