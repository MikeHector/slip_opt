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
imp = perturbColl(dirList);


function [delta_cost] = perturbColl(dirList)
    deltacost = [];
    for i = 1:numel(dirList)
        nameTemp = dirList{i};
        load(nameTemp);
        
        orig_opt_results = opt_results;
        orgSeed = opt_results.X;
        c = opt_results.c;
        
        for k=1:size(orgSeed,1)
            magPerturb = std(orgSeed(k,:));
            orgSeed(k,:) = orgSeed(k,:) + magPerturb * rand(size(orgSeed(k,:)));
        end
        
        perturbSeed = orgSeed; clear orgSeed; clear opt_results;
        
        
        [optimized, opt_results] = RUN_COL(perturbSeed, c, orig_opt_results.apex_velocity,  orig_opt_results.apex_height, 1,  orig_opt_results.end_vel,  orig_opt_results.disturbance_f, NaN);
        
        
        delta_cost{i} =  orig_opt_results.cost - opt_results.cost;
        opt_results.unperturbed = nameTemp;
        opt_results.seed = perturbSeed;
        opt_results.delta_cost = delta_cost{i};
        
        uniqueID = string(datetime, 'dMMyHHmmssSSSS');
        filename = strcat('opt_damping_', uniqueID);
        save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\perturb_results',filename),'opt_results');        
        
    end
end

    
