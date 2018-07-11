% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

delta_damping = 1;
damping_values = 30:delta_damping:400;
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
load('C://Users/mike-/Documents/DRL/collocation/opt_damping_30_baseline.mat')
opt_seed = opt_results.X;

clearvars -except opt_seed apex_vel apex_height damping_values delta_damping

bad_stuff = 0;
too_many_iters = 0;

i = 1;
damping = damping_values(i);
while damping < 300
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel, 0, NaN);
    if (opt_results.flag <= 0) && ((damping_values(i) - damping_values(i - 1)) > 1e-3)
        damping_values(end+1) = (damping_values(i) + damping_values(i-1))/2;
        damping_values = sort(damping_values);           
    else
        filename = strcat('opt_damping_', num2str(cputime*10000000));
        save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\',filename),'opt_results');
%         save(strcat('D:\Documents\DRL\slip_opt\opt_results\damping_results\',filename),'opt_results');
        opt_seed = x_opt_ankle; 
        i = i + 1;
        damping = damping_values(i);
    end
%      damping_values
end

% end
bad_stuff
too_many_iters
