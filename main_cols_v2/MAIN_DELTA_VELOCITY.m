% MHector
% 6.1.18
% Master script to run a bunch of collocations from


%% Increasing
%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

delta_v = .01;
delta_velocity_values = linspace(0, 1.5, (1.5)/delta_v + 1);
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

load('D:\Documents\DRL\slip_opt\opt_results\no_damp_baseline.mat') %My Desktop
% load('C://Users/mike-/Documents/DRL/collocation/opt_damping_30_baseline.mat') %DRL Desktop
% load('C:\Users\Administrator\Documents\DRL\slip_opt\opt_damping_30_baseline.mat') %My laptop
% load('C:\Users\Administrator\Documents\DRL\slip_opt\opt_damping_30_baseline.mat') %DRL laptop
opt_seed = opt_results.X;
damping = opt_results.c;

clearvars -except opt_seed apex_vel apex_height damping delta_velocity_values

bad_stuff = 0;
too_many_iters = 0;

i = 1;
end_velocity = delta_velocity_values(i) + apex_vel;
while end_velocity(i) < 2.6
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, end_velocity, 0, NaN);
%     if (opt_results.flag <= 0) && ((damping_values(i) - damping_values(i - 1)) > 1e-3)
%         damping_values(end+1) = (damping_values(i) + damping_values(i-1))/2;
%         damping_values = sort(damping_values);           
%     else
        filename = strcat('opt_delta_vel_', num2str(cputime*10000000));
%         save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\velocity_results\',filename),'opt_results'); %DRL desktop
        save(strcat('D:\Documents\DRL\slip_opt\opt_results\deltav_results',filename),'opt_results'); %My desktop
%         save(strcat('C:\Users\Administrator\Documents\DRL\slip_opt\opt_results\velocity_results',filename),'opt_results'); %My laptop
%         save(strcat('C:\Users\DRL\Documents\mike_git\slip_opt\opt_results\deltav_results',filename),'opt_results'); %DRL laptop
        opt_seed = x_opt_ankle; 
        i = i + 1;
%         velocity = delta_velocity_values(i);
        end_velocity = delta_velocity_values(i) + apex_vel;
%     end
%      damping_values
end

% end
bad_stuff
too_many_iters

%% Decreasing
%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

delta_v = .01;
delta_velocity_values = linspace(0, -.5, (.5)/delta_v + 1);
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

load('D:\Documents\DRL\slip_opt\opt_results\no_damp_baseline.mat') %My Desktop
% load('C://Users/mike-/Documents/DRL/collocation/opt_damping_30_baseline.mat') %DRL Desktop
% load('C:\Users\Administrator\Documents\DRL\slip_opt\opt_damping_30_baseline.mat') %My laptop
% load('C:\Users\Administrator\Documents\DRL\slip_opt\opt_damping_30_baseline.mat') %DRL laptop
opt_seed = opt_results.X;
damping = opt_results.c;

clearvars -except opt_seed apex_vel apex_height damping delta_velocity_values

bad_stuff = 0;
too_many_iters = 0;

i = 1;
end_velocity = delta_velocity_values(i) + apex_vel;
while end_velocity(i) < 2.6
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, end_velocity, 0, NaN);
%     if (opt_results.flag <= 0) && ((damping_values(i) - damping_values(i - 1)) > 1e-3)
%         damping_values(end+1) = (damping_values(i) + damping_values(i-1))/2;
%         damping_values = sort(damping_values);           
%     else
        filename = strcat('opt_delta_vel_', num2str(cputime*10000000));
%         save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\velocity_results\',filename),'opt_results'); %DRL desktop
        save(strcat('D:\Documents\DRL\slip_opt\opt_results\deltav_results',filename),'opt_results'); %My desktop
%         save(strcat('C:\Users\Administrator\Documents\DRL\slip_opt\opt_results\velocity_results',filename),'opt_results'); %My laptop
%         save(strcat('C:\Users\DRL\Documents\mike_git\slip_opt\opt_results\deltav_results',filename),'opt_results'); %DRL laptop
        opt_seed = x_opt_ankle; 
        i = i + 1;
%         velocity = delta_velocity_values(i);
        end_velocity = delta_velocity_values(i) + apex_vel;
%     end
%      damping_values
end
