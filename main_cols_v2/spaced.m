function [ improve, optKeep ] = spaced(flags, variable)
%spaced returns where improvements are needed or null if no improvements
%are needed
%   The optimization must have N variables with positive
%   flags that are +/- M % away from being evenly spaced
    N = 20;
    M = 1; %percent
    middle = floor(mean(variable));
    %Find max variable
    badCounter = 0;
    maxVar = [];
    for i = middle:length(variable)
        if flags(i) < 0
            badCounter = badCounter + 1;
        end
        if badCounter == 5
            maxVar = variable(i - 5);
            break
        elseif badCounter < 5 && flags(i) > 0
            badCounter = 0;
        end
    end
    if maxVar == []
        maxVar = variable(i);
    end        
    
    
    %Find min variable
    badCounter = 0;
    minVar = nan;
    for i = linspace(middle, 1, middle)
        if flags(i) < 0
            badCounter = badCounter + 1;
        end
        if badCounter == 5
            minVar = variable(i + 5);
            break
        elseif badCounter < 5 && flags(i) > 0
            badCounter = 0;
        end
    end
    if isnan(minVar)
        minVar = variable(i);
    end
    
    dist = linspace(minVar, maxVar,N);
    span = abs(minVar - maxVar);
    improve = [];
    optKeep = [];
    for i = 1:length(dist)
        highBound = dist(i) + M/100 * span;
        lowBound = dist(i) - M/100 * span;
        [~, highSet] = find(variable < highBound);
        [~, lowSet] = find(variable > lowBound);
        goodSet = intersect(highSet, lowSet);
        goodFlags = flags(goodSet);
        if numel(goodSet) > 0 && max(goodFlags) > 0
            %We gucci
            optKeep = [optKeep; max(goodSet)];
        else
            improve = [improve; highBound, lowBound];
        end
    end

end