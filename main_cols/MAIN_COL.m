% MHector
% 7/25/18
% Master script to run a bunch of collocations from

%Let's run over several variable values using previous variable value
%solution as a seed for the next optimization
clear; clc;

%for varName = {'c', 'deltav', 'FDisturb', 'TDA', 'apexVel'}
% colStruc = getColStruc('varName');

colStruc.direction = {'up','down'};
colStruc.varName = 'c';
colStruc.deltaVar = 1;
colStruc.varMax = 500;
colStruc.varMin = 0;

for k = 1:numel(colStruc.direction)
    direction = colStruc.direction{k};
    
    %Load baseline
    load('THE BASELINE WILL GO HERE')
    opt_seed = opt_results.X;
    param = opt_results.param;
    clear opt_results

    %Define variables to be changed during analysis
    varName = 'c';
    deltaVar= 1;
    var = 0;

    infeasibleCounter = 0;

    while infeasibleCounter < 20 && var <= colStruc.maxVar && var >= colStruc.minVar
        [new_opt_seed, opt_results] = RUN_COL(opt_seed, param);

        %Save the coll
        uniqueID = string(datetime, 'dMMyHHmmssSSSS');
        filename = strcat('opt_damping_', varName, '_', uniqueID);
        save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\',filename),'opt_results');
        opt_seed = new_opt_seed; 

        %Track stuff
        if opt_results.flag < 0
            infeasibleCounter = infeasibleCounter + 1;
        end
        
        %Increment variable
        param = opt_results.param;
        if strcmp(direction,'up')
            param.(var) = opt_results.param.(var) + deltaVar;
        elseif strcmp(direction, 'down')
            param.(var) = opt_results.param.(var) - deltaVar;
        end
        
        %Clean
        clear new_opt_seed opt_results
    end
end
