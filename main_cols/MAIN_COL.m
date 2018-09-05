% MHector
% 9/4/18
% Master script to run a bunch of collocations from

%Let's run over several variable values using previous variable value
%solution as a seed for the next optimization
clear; clc;
colStrucArray = ColStrucBuilder();
% colStrucArray = ColStrucBuilderTest();
fieldNames = fieldnames(colStrucArray);

saveDir = getSaveDir('Michael-PC');

for m = [8]
    colStruc = colStrucArray.(fieldNames{m});

    for k = 1:numel(colStruc.direction)
        direction = colStruc.direction{k};

        %Load baseline seed and parameters
        load('baseline4R')
%         load('C:\Users\mike-\Documents\DRL\slip_opt\opt_results\opt_disturbance_f_180820182234015610.mat')
%         load( 'C:\Users\mike-\Documents\DRL\collocation\opt_results\opt_R_leg_150820181810270730.mat')
        opt_seed = opt_results.X;
        param = opt_results.param;
        param.(colStruc.varName) = colStruc.var;
        clear opt_results

        NumPerturb = 5;
        badCounterGlobe = 0;
        tryCounter = 1;
        opt_results.param.flag = 0;
        while badCounterGlobe < 7 && param.(colStruc.varName) <= colStruc.varMax && param.(colStruc.varName) >= colStruc.varMin
            opt_results.param.flag = 0;
            while tryCounter <= NumPerturb && opt_results.param.flag <= 0
                [DV_out, opt_results] = RUN_COL2(opt_seed, param, 0); %softmax
                
                %Perturb if it hasn't converged
                if opt_results.param.flag <= 0 && tryCounter < NumPerturb
                    for i = 1:size(opt_results.X,1)
                         opt_seed(i,:) = opt_results.X(i,:) +...
                         .5 * std(opt_results.X(i,:)) * rand(size(opt_seed(2,:)));
                    end
                end
                tryCounter = tryCounter + 1;
            end
            tryCounter = 0;
                

            %Save the coll
            uniqueID = string(datetime, 'dMMyHHmmssSSSS');
            filename = strcat('opt_', colStruc.varName, '_', uniqueID);
            save(strcat(saveDir,filename),'opt_results');

            %Save optimized decision variables as new seed
            opt_seed = DV_out; 

            %Track infeasible counter
            if (opt_results.param.flag <= 0)
                badCounterGlobe = badCounterGlobe + 1;
            elseif opt_results.param.flag > 0
                badCounterGlobe = 0; %Reset counter if not infeasible/ bad
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