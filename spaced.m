function [ output_args ] = spaced(flags, variable)
%spaced returns where improvements are needed or null if no improvements
%are needed
%   The optimization must have N variables with positive
%   flags that are +/- M % away from being evenly spaced
    N = 20;
    M = 1; %percent
    dist = linspace(min(variable), max(variable),N);
    span = abs(min(variable) - max(variable));
    improve = [];
    for i = 1:length(dist)
        highBound = dist(i) + M/100 * span;
        lowBound = dist(i) - M/100 * span;
        if (any(variable > lowBound) && any(variable < highBound)) && flags(i) > 0
            %We gucci
        else
            improve = [improve; highBound, lowBound];
        end
    end

end