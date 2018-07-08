function [ improve ] = spaced(flags, variable)
%spaced returns where improvements are needed or null if no improvements
%are needed
%   The optimization must have N variables with positive
%   flags that are +/- M % away from being evenly spaced
    N = 20;
    M = 1; %percent
    middle = mean(variable);
    %Find max variable
    badCounter = 1;
    for i = middle:length(variable)
        if flag(i) < 0
            badCounter = badCounter + 1;
        end
        if badCounter == 5
            maxVar = variable(i -5);
            break
        end
    end
    
    %Find min variable
    badCounter = 1;
    for i = 1:middle
        if flag(i) < 0
            badCounter = badCounter + 1;
        end
        if badCounter == 5
            maxVar = variable(i -5);
            break
        end
    end
    
    dist = linspace(minVar, maxVar,N);
    span = abs(minVar - maxVar);
    improve = [];
    for i = 1:length(dist)
        highBound = dist(i) + M/100 * span;
        lowBound = dist(i) - M/100 * span;
        [~, highSet] = find(variable < highBound);
        [~, lowSet] = find(variable > lowBound);
        goodSet = intersect(highSet, lowSet);
        goodFlags = flags(goodSet);
        if numel(goodSet) > 0 && max(goodFlags) > 0
            %We gucci
        else
            improve = [improve; highBound, lowBound];
        end
    end

end