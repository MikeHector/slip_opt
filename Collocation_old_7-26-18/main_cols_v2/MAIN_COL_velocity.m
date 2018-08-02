% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

%% Increasing
delta_vel = .02;
vel_values = linspace(1, .5, (1-.5)/delta_vel);
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

load('D:\Documents\DRL\slip_opt\main_cols_v2\opt_damping_150.mat')
opt_seed = opt_results.X;
damping = opt_results.c;

clearvars -except opt_seed apex_vel apex_height damping vel_values delta_vel

bad_stuff = 0;
too_many_iters = 0;

i = 1;
for i =1:length(vel_values)
    ankles_on = 1;
    apex_vel = vel_values(i);
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel, 0, atan2(opt_seed(2,1),opt_seed(1,1)));
%     if opt_results.flag <= 0
%         vel_values(end+1) = (vel_values(i) + vel_values(i-1))/2;
%         vel_values = sort(vel_values);
%     else
    filename = strcat('opt_vel_', num2str(cputime*10000000));
    save(strcat('D:\Documents\DRL\slip_opt\opt_results\velocity_results\',filename),'opt_results');
    opt_seed = x_opt_ankle; 
%     end

     
%      if abs(min(diff(vel_values))) < 1e-6
%          break
%      end
%      vel_values
     if opt_results.flag == 0
         too_many_iters = too_many_iters + 1;
     elseif opt_results.flag < 0
         bad_stuff = bad_stuff + 1;
     end
end

% end
bad_stuff
too_many_iters

%% Decreasing
delta_vel = .02;
vel_values = linspace(1, 3, 2/delta_vel);
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

load('D:\Documents\DRL\slip_opt\main_cols_v2\opt_damping_150.mat')
opt_seed = opt_results.X;
damping = opt_results.c;

clearvars -except opt_seed apex_vel apex_height damping vel_values delta_vel

bad_stuff = 0;
too_many_iters = 0;


for i =1:length(vel_values)
    ankles_on = 1;
    apex_vel = vel_values(i);
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel, 0, atan2(opt_seed(2,1),opt_seed(1,1)));
%     if opt_results.flag <= 0
%         vel_values(end+1) = (vel_values(i) + vel_values(i-1))/2;
%         vel_values = sort(vel_values);
%     else
    filename = strcat('opt_vel_', num2str(cputime*10000000));
    save(strcat('D:\Documents\DRL\slip_opt\opt_results\velocity_results\',filename),'opt_results');
    opt_seed = x_opt_ankle; 
%     end

     
%      if abs(min(diff(vel_values))) < 1e-6
%          break
%      end
%      vel_values
     if opt_results.flag == 0
         too_many_iters = too_many_iters + 1;
     elseif opt_results.flag < 0
         bad_stuff = bad_stuff + 1;
     end
end

% end
bad_stuff
too_many_iters