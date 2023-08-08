function const = constants()
% CONSTANTS - returns a structure containing constants used in the
% across the scripts

const.areanames = ["HPC", "PFC"];

% DEPRECATED : use the Option struct's .shortcut field
const.HPC = 1;
const.PFC = 2;
const.THETA = 1;
const.DELTA = 2;
const.RIPPLE = 3;

const.directions = ["hpc-hpc", "hpc-pfc"];

% NOTE: double check which used for prev papes
% animals = ["JS21","ZT2","ER1","JS14","JS13","JS17"];
const.all_animals = ["ER1", "JS13", "JS14", "JS15", "JS17", "JS21", "ZT2"];

% NOTE: USE BETTER COLORS
const.hpccolor = 'blue';
[~,const.pfccolor] = colornames('Wikipedia', 'orange');
