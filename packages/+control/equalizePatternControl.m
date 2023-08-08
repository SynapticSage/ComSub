function [equalized, warned] = equalizePatternControl(cellOfWindows, Option)
% equalized the number of windows between a pattern and its control
% cellOfWindows has to be formatted so that the first nPatterns are the
% patterns and the second half are the controls

equalized = cellOfWindows;
warned = false;

for pattern = 1:Option.nPatterns

    mods = mod(1-pattern:length(cellOfWindows)-pattern, Option.nPatterns);
    mods = find(mods == 0);
    mods = setdiff(mods, pattern);
    for m = mods

        disp("equalizing pattern " + pattern + " and control " + m);

        curr_pattern = cell2mat(cellOfWindows(:,pattern));
        curr_control = cell2mat(cellOfWindows(:,m));
        
        % if control window pattern is empty (due to removing overlaps), use
        % the original number of the windows in the pattern.
        if isempty(curr_control)
            warning("control" + pattern +...
                    " is empty,using original number of pattern windows");
            warned = true;
            break;
        end
        
        [real_pattern, real_control] = ...
            windows.removeUntilEqual(curr_pattern, curr_control, "matched", 50);
        equalized{pattern} = real_pattern;
        equalized{m}       = real_control;
    end
end


end

