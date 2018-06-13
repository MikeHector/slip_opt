% MHector
% 6.1.18
% Master script to run a bunch of collocations from
%Something new
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
load('opt_f_disturb_baseline') %Baseline from other analysis
seed = opt_results.X;
apex_vel = opt_results.apex_velocity;
end_vel = opt_results.apex_velocity;
apex_height = opt_results.apex_height;
damping_values = opt_results.c;
disturbance_f = 0;
baseline_TD = atan2(opt_results.y(1), opt_results.x(1));
TD_disturb = linspace(0,-.2,15);
TD_angle = baseline_TD + TD_disturb;
clear opt_results

bad_stuff = 0;
too_many_iters = 0;
wb = waitbar(0, 'Optimizing unoptimized optimizations');
%Start from seed and go up
for i = 1:length(TD_angle)
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel, apex_height, ankles_on, end_vel, disturbance_f, TD_angle(i));
    save(strcat('opt_TD_disturb_', num2str(cputime*10000000)),'opt_results');
    seed = x_opt_ankle;
    if opt_results.flag < 0
        bad_stuff = bad_stuff + 1;
    elseif opt_results.flag == 0
        too_many_iters = too_many_iters + 1;
    end
    close(wb)
    wb = waitbar(i/length(TD_angle),'Optimizing unoptimized optimizations');
end
% figure; plot(opt_results.x, opt_results.y)
% figure; plot(opt_results.t, opt_results.Tleg, 'b'); hold on; plot(opt_results.t, opt_results.Tankle, 'r')

bad_stuff
too_many_iters
