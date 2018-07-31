% MHector
% 7/25/18
% Master script to run a bunch of collocations from

%Let's run over several variable values using previous variable value
%solution as a seed for the next optimization
clear; clc;
colStrucArray = ColStrucBuilder();
% colStrucArray = ColStrucBuilderTest();
fieldNames = fieldnames(colStrucArray);

for m = 2:numel(fieldNames)
    colStruc = colStrucArray.(fieldNames{m});
    

    for k = 1:numel(colStruc.direction)
        direction = colStruc.direction{k};

        %Load baseline seed and parameters
        load('opt_c_300720181641270440')
        opt_seed = opt_results.X;
        param = opt_results.param;
        param.(colStruc.varName) = colStruc.var;
        clear opt_results

        infeasibleCounter = 0;

        while infeasibleCounter < 20 && param.(colStruc.varName) <= colStruc.varMax && param.(colStruc.varName) >= colStruc.varMin
            [DV_out, opt_results] = RUN_COL(opt_seed, param);

            %Save the coll
            uniqueID = string(datetime, 'dMMyHHmmssSSSS');
            filename = strcat('opt_', colStruc.varName, '_', uniqueID);
            save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\',filename),'opt_results');

            %Save optimized decision variables as new seed
            opt_seed = DV_out; 

            %Track infeasible counter
            if opt_results.param.flag < 0
                infeasibleCounter = infeasibleCounter + 1;
            elseif opt_results.param.flag >= 0
                infeasibleCounter = 0; %Reset counter if not infeasible
            end

            %Increment collocation variable
            param = opt_results.param;
            if strcmp(direction,'up')
                param.(colStruc.varName) = opt_results.param.(colStruc.varName) + colStruc.deltaVar;
            elseif strcmp(direction, 'down')
                param.(colStruc.varName) = opt_results.param.(colStruc.varName) - colStruc.deltaVar;
            end

            %Clean
            clear DV_out opt_results
        end
    end
end