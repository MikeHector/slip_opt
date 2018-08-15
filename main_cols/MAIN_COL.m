% MHector
% 7/25/18
% Master script to run a bunch of collocations from

%Let's run over several variable values using previous variable value
%solution as a seed for the next optimization
clear; clc;
colStrucArray = ColStrucBuilder();
% colStrucArray = ColStrucBuilderTest();
fieldNames = fieldnames(colStrucArray);

saveDir = getSaveDir('DRL-PC');

for m = [8]
    colStruc = colStrucArray.(fieldNames{m});

    for k = 1:numel(colStruc.direction)
        direction = colStruc.direction{k};

        %Load baseline seed and parameters
%         load('baseline4')
        load('baseline_no_electrical_loses.mat')
        opt_seed = opt_results.X;
        param = opt_results.param;
        param.(colStruc.varName) = colStruc.var;
        clear opt_results

        badCounter = 0;

        while badCounter < 20 && param.(colStruc.varName) <= colStruc.varMax && param.(colStruc.varName) >= colStruc.varMin
            [DV_out, opt_results] = RUN_COL(opt_seed, param);

            %Save the coll
            uniqueID = string(datetime, 'dMMyHHmmssSSSS');
            filename = strcat('opt_', colStruc.varName, '_', uniqueID);
            save(strcat(saveDir,filename),'opt_results');

            %Save optimized decision variables as new seed
            opt_seed = DV_out; 

            %Track infeasible counter
            if (opt_results.param.flag <= 0)
                badCounter = badCounter + 1;
            elseif opt_results.param.flag > 0
                badCounter = 0; %Reset counter if not infeasible/ bad
            end

            %Increment collocation variable
            param = opt_results.param;
            if strcmp(direction,'up')
                param.(colStruc.varName) = opt_results.param.(colStruc.varName) + colStruc.deltaVar;
            elseif strcmp(direction, 'down')
                param.(colStruc.varName) = opt_results.param.(colStruc.varName) - colStruc.deltaVar;
            end
            disp([colStruc.varName , ' has been incremented to ',num2str(param.(colStruc.varName))]);
            %Clean
            clear DV_out opt_results
        end
    end
end