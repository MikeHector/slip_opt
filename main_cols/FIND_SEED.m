% MHector
% 7/25/18
% Master script to run a bunch of collocations from

%Let's run over several variable values using previous variable value
%solution as a seed for the next optimization
clear; clc;
% colStrucArray = ColStrucBuilder();
% colStrucArray = ColStrucBuilderTest();
% fieldNames = fieldnames(colStrucArray);

newFileName = 'baseline5';
%Load baseline seed and parameters
load('baseline5')
opt_seed = opt_results.X;
param = opt_results.param;

iterationCounter = 0;
lowest_cost = 1e6;

while iterationCounter < 20000
    [DV_out, opt_results] = RUN_COL(opt_seed, param);

    if opt_results.cost < lowest_cost && opt_results.param.flag > 0
        lowest_cost = opt_results.cost;
        %Save the coll
%         save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\',filename),'opt_results');
        save(strcat('C:\Users\mike-\Documents\DRL\collocation\main_cols\',newFileName),'opt_results');
        %Perturb it a little
        for i = 1:size(opt_results.X,1)
            opt_seed(i,:) = opt_results.X(i,:) +...
                .5 *std(opt_results.X(i,:)) * rand(size(opt_seed(2,:)));
        end
    else
        %Load best seed and
        %Perturb it a little more
        load(strcat('C:\Users\mike-\Documents\DRL\collocation\main_cols\',newFileName))
        for i = 1:size(opt_results.X,1)
            opt_seed(i,:) = opt_results.X(i,:) +...
                .5 * std(opt_results.X(i,:)) * rand(size(opt_seed(2,:)));
        end
    end

    %Increment counter
    iterationCounter = iterationCounter + 1;

    %Clean
    clear DV_out opt_results
end