function nonOverlappingpattern = removeOverlapBetweenPattern...
                                  (windowsToPreserve, windowsToRemove)
    %
    % This function removes windows from windowsToRemove that overalp windows
    % in windowsToPreserve. This in essence, can remove theta control windows
    % who have times in theta windows.
    %
    % patternToPreserve is a cell of start&stop times for a single pattern,
    % e.g. theta.
    % returns a cell with start&stop times non-overlapping
    
    % for i = 1:nPatterns
    %     curr_controlPattern = controlWindows{i};
    %     curr_Pattern = cellOfWindows{i};
    %     curr_ControlPattern = curr_controlPattern(1:1000,:);
    %     curr_Pattern = curr_Pattern(1:200,:);
    %     test_preserve{i} = curr_Pattern;
    %     test_remove{i} = curr_ControlPattern;
    % end
    %
    tic
    [nPossibleWindows,~] = size(windowsToRemove);
    [nExistingWindows,~] = size(windowsToPreserve);
    curr_existingwindow = zeros(1,2); % why without this it does not work??
    for j = 1:nExistingWindows
        curr_existingwindow = windowsToPreserve(j,:);
        tracker = 1;
        while tracker <= nPossibleWindows
            curr_Possiblewindow = windowsToRemove(tracker,:);
            if windows.detectOverlap (curr_existingwindow(1), curr_existingwindow(2), ...
                    curr_Possiblewindow(1),curr_Possiblewindow(2))
                windowsToRemove(tracker,:) = [];
                [nPossibleWindows,~] = size(windowsToRemove);
            else
                tracker = tracker+1;
            end
        end
    end
    nonOverlappingpattern = windowsToRemove;
    disp("removeOverlapBetweenPattern took: " + toc + " seconds");
end
