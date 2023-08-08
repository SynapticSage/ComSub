function [new_pattern1, new_pattern2]=removeUntilEqual(pattern1, pattern2, method, binsToMatch)

% this function takes in two pattern windows and remove the one with
% more windows until they are equal.

% if method is "random", it randomly remove windows from the pattern
% with more windows
% if the method is "matched", it will

[n1,~] = size(pattern1);
[n2,~] = size(pattern2);

if n1>n2
    toRemove = pattern1;
else
    toRemove = pattern2;
end

if method == "random"
    removeSize = abs(n1-n2);
    
    % indices = randperm(removeSize); 2023 remove
    indices = randperm(max(n1, n2), removeSize); % 2023
    toRemove(indices',:) = [];
    indices = sort(indices);
    
    if n1>n2
        new_pattern1 = toRemove;
        new_pattern2 = pattern2;
    else
        new_pattern2 = toRemove;
        new_pattern1 = pattern1;
    end
    
elseif method == "matched"
    timeaxis = [min(pattern1(1,1), pattern2(1,1)),max(pattern1(end,2),pattern2(end,2))];
    timeBins = linspace(timeaxis(1), timeaxis(2),binsToMatch+1);
    
    p1Tracker = 1;
    p2Tracker = 1;
    
    p1ByBins = cell(1,binsToMatch);
    p2ByBins = cell(1,binsToMatch);
    
    % select windows of the two patterns that belong to each time bins
    for i = 1:binsToMatch
        startTime = timeBins(i);
        endTime = timeBins(i+1);
        
        for j = p1Tracker+1:n1
            if pattern1(j,1)>endTime
                p1ByBins{i} = pattern1(p1Tracker+1:j, :);
                p1Tracker = j;
                break
            end
        end
        
        for j = p2Tracker+1:n2
            if pattern2(j,1)>endTime
                p2ByBins{i} = pattern2(p2Tracker+1:j, :);
                p2Tracker = j;
                break
            end
        end
%         disp(p2Tracker)
    end
    
    % for each of the bins, select the minimum of the two patterns, and
    % randomly remove the others until they are equal in number
    
    for i = 1:binsToMatch
        numP1Windows = size(p1ByBins{i},1);
        numP2Windows = size(p2ByBins{i},1);
        
        if numP1Windows > numP2Windows
            if numP2Windows >0
                randIdcs = randperm(numP1Windows,numP2Windows);
                toExclude = setdiff(1:numP1Windows, randIdcs);
%                 assert(size(p1ByBins{i}(~randIdcs, :),1 )==(numP1Windows - numP2Windows))
                p1ByBins{i}(toExclude, :) = [];
            else
                p1ByBins{i} = [];
            end
        elseif numP1Windows < numP2Windows
            if numP1Windows > 0
                randIdcs = randperm(numP2Windows, numP1Windows);
                toExclude = setdiff(1:numP2Windows, randIdcs);
%                 assert(numel(randIdcs) == numP1Windows)
                p2ByBins{i}(toExclude,  :) = [];
            else
                p2ByBins{i} = [];
            end
        end
        
    end
    % finally, concatenate all the remaining windows into a windows matrix
    new_pattern1 = [];
    new_pattern2 = [];
    for i = 1:binsToMatch
        new_pattern1 = [new_pattern1; p1ByBins{i}];
        new_pattern2 = [new_pattern2; p2ByBins{i}];
    end
end
end
