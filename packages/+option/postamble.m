function Option = postable(Option)

%% ===============================================
%% Shortcut/alias variables to improve readability
%% ===============================================
animal = Option.animal;
Option.shortcut.THETA  = 1;
Option.shortcut.DELTA  = 2;
Option.shortcut.RIPPLE = 3;
if Option.sourceArea == "CA1"
    Option.shortcut.HPC = 1;
    Option.shortcut.PFC = 2;
else
    Option.shortcut.PFC = 1;
    Option.shortcut.HPC = 2;
end
Option.patternNames = ["theta", "delta", "ripple"];
winSize = Option.winSize(1);
Option.frequenciesPerPattern = [6 10; 0.5 4; 150 200];
[Option.nPatterns,~] = size(Option.frequenciesPerPattern);
if Option.singleControl == true
    Option.nPatternAndControl = Option.nPatterns+1;
else
    Option.nPatternAndControl = Option.nPatterns*2;
end

if Option.midpattern
    Option.nPatternAndControl = Option.nPatternAndControl + Option.nPatterns;
end

if ~any(contains(Option.patternNames,"control"))
    Option.patternNames = [Option.patternNames; Option.patternNames+"-control"]';
    Option.patternNames = Option.patternNames(:)';
end
if Option.midpattern
    Option.patternNames = [Option.patternNames; Option.patternNames+"-mid"]';
    Option.patternNames = Option.patternNames(:)';
    figuredefine("-clearpermfolder")
    figuredefine("-permfolder", "midpattern=true")
    disp("Set figuredefine to " + figuredefine("-getpermfolder"))
else
    figuredefine("-clearpermfolder")
    disp("Set figuredefine to " + figuredefine("-getpermfolder"))
end

if any(contains(Option.generateH, "fromWpli"))
    Option.patternNamesFull = Option.patternNames + "-wpli";
    Option.genH_name = "wpli";
elseif any(contains(Option.generateH, "fromCoherence"))
    Option.patternNamesFull = Option.patternNames + "-coherence";
    Option.genH_name = "coherence";
elseif any(contains(Option.generateH, ["fromSpectra", "fromHilbert"]))
    Option.patternNamesFull = Option.patternNames + "-power";
    Option.genH_name = "power";
elseif any(contains(Option.generateH, "fromCP"))
    Option.patternNamesFull = Option.patternNames + "-cp";
    Option.genH_name = "c_given_p";
else
    error("No valid Option.generateH")
end


% if Option.genH_name == "power"
%     Option.nPatternAndControl = Option.nPatternAndControl - 1;
%     Option.patternNamesFull   = Option.patternNamesFull(1:end-1);
%     Option.patternNames       = Option.patternNames(1:end-1);
% end
