function [saveDir] = getSaveDir(TypeZeroForListOfInputs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if TypeZeroForListOfInputs == 0
        disp('Place -> Input')
        disp('My DRL computer -> DRL-PC')
        disp('My Home computer -> Alluring-PC')
        disp('Michaels DRL computer -> Michael-PC')

    elseif strcmp(TypeZeroForListOfInputs, 'DRL-PC')
        saveDir = 'C:\Users\mike-\Documents\DRL\collocation\opt_results\';

    elseif strcmp(TypeZeroForListOfInputs, 'Alluring-PC')
        saveDir = ' ';

    elseif strcmp(TypeZeroForListOfInputs, 'Michael-PC')
        saveDir = ' ';
    
    end



end

