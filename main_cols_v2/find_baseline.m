% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

delta_damping = 1;
damping_values = 0;
% damping_values = [damping_values, 200:2:1500];
apex_vel = 1; apex_height = 1.1; 

% %Null seed
% x = linspace(-.15, .15, 40);
% y = linspace(0, 0, 40);
% r0 = linspace(.9, .9, 40);
% dx = linspace(.4,.4,40);
% dy = linspace(0,0,40);
% dr0 = linspace(0,0,40);
% Tl = linspace(15,15,40);
% Ta = linspace(0,0,40);
% t = linspace(0,.2,40);
% seedy = [x; y; r0; dx; dy; dr0; Tl; Ta; t];

% load('D:\Documents\DRL\slip_opt\opt_results\no_damp_baseline.mat') 
% load('C://Users/mike-/Documents/DRL/collocation/opt_results/baselines/no_damp_baseline.mat')
load('C:\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\new_obj_func\opt_damping_baseline_60.mat')
% load('C://Users/mike-/Documents/DRL/collocation/opt_damping_30_baseline.mat')
% load('C:\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_318924375000.mat')
% load('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_363118125000.mat')
opt_seed = opt_results.X;

clearvars -except opt_seed apex_vel apex_height damping_values delta_damping

bad_stuff = 0;
too_many_iters = 0;

i = 1;
damping = damping_values(i);
bad_counter = 0;
count = 1;
lowest_cost = 1e6;
save_counter = 0;
while count < 50
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel, 0, NaN);
%     if (opt_results.flag <= 0) && ((damping_values(i) - damping_values(i - 1)) > 1e-3)
%         damping_values(end+1) = (damping_values(i) + damping_values(i-1))/2;
%         damping_values = sort(damping_values);           
%     else
%         uniqueID = string(datetime, 'dMMyHHmmssSSSS');
        cost_track(count) = opt_results.cost;
        if opt_results.flag > 0 && opt_results.cost < lowest_cost && save_counter <= 20 && opt_results.cost > 0
            lowest_cost = opt_results.cost;
            uniqueID = 'baseline_Rs';
            filename = strcat('opt_damping_', uniqueID);
            save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\new_obj_func\',filename),'opt_results');
%         save(strcat('D:\Documents\DRL\slip_opt\opt_results\damping_results\',filename),'opt_results');
            opt_seed = x_opt_ankle;
            
            %Perturb it a little
            for k=1:size(opt_seed,1)
                magPerturb = .1 * std(opt_seed(k,:));
                opt_seed(k,:) = opt_seed(k,:) + magPerturb * rand(size(opt_seed(k,:)));
            end
            save_counter = save_counter + 1;
        
        else 
            %Load best so far
            load('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\new_obj_func\opt_damping_baseline.mat')
            opt_seed = opt_results.X;
            %Perturb it more
            for k=1:size(opt_seed,1)
                magPerturb = .5 * std(opt_seed(k,:));
                opt_seed(k,:) = opt_seed(k,:) + magPerturb * rand(size(opt_seed(k,:)));
            end
            save_counter = 0;
        end

    count = count + 1
       
end

% end
bad_stuff
too_many_iters
