% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

% Load initial seed
load('opt_damping_147') %Baseline some damping
seed = opt_results.X;
apex_vel = opt_results.apex_velocity;
damping_values = opt_results.c;
apex_height = opt_results.apex_height;
bad_stuff = 0;
too_many_iters = 0;
deltav = linspace(0,1,30);

%Start from seed and go down a bit
for i = 1:length(deltav)
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel, apex_height, ankles_on, apex_vel + deltav(i));
    save(strcat('opt_deltav_', num2str(cputime*10000000)),'opt_results');
    seed = x_opt_ankle; 
    if opt_results.flag < 0
        bad_stuff = bad_stuff + 1;
    elseif opt_results.flag == 0
        too_many_iters = too_many_iters + 1;
    end
end

%Load initial seed
clear
load('opt_damping_147') %Baseline some damping
seed = opt_results.X;
apex_vel = opt_results.apex_velocity;
apex_height = opt_results.apex_height;
damping_values = opt_results.c;
clear opt_results
deltav = linspace(0,-.5,15);
bad_stuff = 0;
too_many_iters = 0;

for i = 1:length(deltav)
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel, apex_height, ankles_on, apex_vel + deltav(i));
    save(strcat('opt_deltav_', num2str(cputime*10000000)),'opt_results');
    seed = x_opt_ankle; 
    if opt_results.flag < 0
        bad_stuff = bad_stuff + 1;
    elseif opt_results.flag == 0
        too_many_iters = too_many_iters + 1;
    end
end

bad_stuff
too_many_iters
