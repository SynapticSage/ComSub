function equalized = equalizeWindowsAcrossPatterns(cellOfWindows)
% this function equalize the windows across patterns within a cell of
% windows
  
    least_windows = cellOfWindows{1};
    location = 1;
    for i = 1:length(cellOfWindows)
        if size(cellOfWindows{i},1) < size(least_windows,1)
            least_windows = cellOfWindows{i};
            location = i;
        end
    end
    
    for i = 1:length(cellOfWindows)
        if i == location
            continue;
        end
        [~, equalToLeast] = windows.removeUntilEqual(least_windows, ...
                                 cellOfWindows{i},"matched", 50);
        cellOfWindows{i} = equalToLeast;
    end
    
    equalized = cellOfWindows;
end

