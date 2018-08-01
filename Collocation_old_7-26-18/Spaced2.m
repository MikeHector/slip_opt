%M Hector
%7.17.18
%Spaced2
% clear; clc;
% directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt*';
% variable_name = 'c';
% n = 50;
% m = 1; %percent
% [g i] = space_func(directory, variable_name, n, m)

function [golden, improve] = space_func(directory,variable_name, n, m)
    strucc = dir(directory); 
    for i = 1:length(strucc)
        filename = strucc(i).name;
        load(strcat(directory(1:end-4), filename))
        res.var = opt_results.(variable_name); res.flag = opt_results.flag;
        res.filename = filename;
        S{i} = res;
        if res.flag >= 0 %Make sure they aren't infeasible
            var_for_max(i) = res.var;
        end
            
    end
    maxVar = max(var_for_max);
    minVar = min(var_for_max);
    mVar = (maxVar - minVar) * .01;
    nArray = linspace(minVar,maxVar,n);
    
    for q = 1:numel(S)
        fuckme(q) = S{q}.var;
    end
    
    improve_count = 1;
    j = 1;
    for i = 1:numel(nArray)
        found_one = false;

        k = 1;
        while ~found_one && k <= numel(S)
            if S{k}.var < nArray(i) + mVar && S{k}.var > nArray(i) - mVar && S{k}.flag > 0
                golden{j} = S{k};
                j = j + 1;
                found_one = true;
            else
                k = k + 1;
%                 if k > numel(S)
%                     k = k - 1;
%                 end
            end
        end
        
        if k == numel(S) + 1
            k = k - 1;
        end
        
        if ~found_one
            improve(improve_count) = nArray(i);
            improve_count = improve_count + 1;
        end
    end

end
