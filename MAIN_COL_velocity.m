% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;


% apex_vel = 1; apex_height = 1.1; 

%Load initial seed
% load('opt_damping_147') %Baseline some damping

% bad_stuff = 0;
% too_many_iters = 0;
% apex_vel_values = linspace(apex_vel,.25,30);
% %Start from seed and go down a bit
% for i = 1:length(apex_vel_values)
%     ankles_on = 1;
%     [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel_values(i), apex_height, ankles_on);
%     save(strcat('opt_velocity_', num2str(cputime*10000000)),'opt_results');
%     seed = x_opt_ankle; 
%     if opt_results.flag < 0
%         bad_stuff = bad_stuff + 1;
%     elseif opt_results.flag == 0
%         too_many_iters = too_many_iters + 1;
%     end
% endseed = opt_results.X;
% apex_vel = opt_results.apex_velocity;
% damping_values = opt_results.c;
% clear opt_results


%Load initial seed
load('opt_velocity_701494843750') %Baseline from other analysis
seed = opt_results.X;
apex_vel = opt_results.apex_velocity;
apex_height = opt_results.apex_height;
damping_values = opt_results.c;
clear opt_results
apex_vel_values = linspace(apex_vel,6,15);
bad_stuff = 0;
too_many_iters = 0;
%Start from seed and go up
for i = 1:length(apex_vel_values)
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel_values(i), apex_height, ankles_on);
    save(strcat('opt_velocity_', num2str(cputime*10000000)),'opt_results');
    seed = x_opt_ankle; 
    if opt_results.flag < 0
        bad_stuff = bad_stuff + 1;
    elseif opt_results.flag == 0
        too_many_iters = too_many_iters + 1;
    end
end

bad_stuff
too_many_iters
