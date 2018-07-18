%MHector
%7.17.18
%Pertubations
clc; clear;
parentDir = 'C:\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt*';
[good, ~] = space_func(parentDir, 'c', 50, 1);
dirList = [];
for i= 1:numel(good)
    dirList{i} = strcat(parentDir(1:end-4), good{i}.filename);
end


function [improvements] = perturbColl(dirList)
    for i = 1:numel(dirList)
        nameTemp = dirList{i}
    end
end

    
