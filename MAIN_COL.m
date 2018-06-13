% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

damping_values = 0:1:199;
damping_values = [damping_values, 200:2:1500];
apex_vel = 1; apex_height = 1.1; 

%Load initial seed
load('opt_030') %Baseline no damping which has been check from multiple seeds
seed = optimized; 
clearvars -except seed apex_vel apex_height damping_values

bad_stuff = 0;
too_many_iters = 0;
for i = 1:length(damping_values)
    
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values(i), apex_vel, apex_height, ankles_on);
    save(strcat('opt_damping_', num2str(damping_values(i))),'opt_results');
    seed = x_opt_ankle; 
    if opt_results.flag < 0
        bad_stuff = bad_stuff + 1;
    elseif opt_results.flag == 0
        too_many_iters = too_many_iters + 1;
    end

end
bad_stuff
too_many_iters
