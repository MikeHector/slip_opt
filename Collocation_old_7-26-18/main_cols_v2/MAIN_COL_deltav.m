% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

dv = .05;
vel_values = 0:dv:2;
% damping_values = [damping_values, 200:2:1500];
apex_vel = 1; apex_height = 1.1; 
damping = 40;

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

load('C://Users/mike-/Documents/DRL/collocation/opt_results/damp_baseline.mat')
opt_seed = optimized;

clearvars -except opt_seed apex_vel apex_height damping_values delta_damping

bad_stuff = 0;
too_many_iters = 0;

i = 1;
% apex_vel = vel_values(i);
while vel_values(i) <= 2
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel + vel_values(i), 0, atan2(opt_seed(2,1),opt_seed(1,1)));
    if opt_results.flag <= 0
        vel_values(end+1) = (vel_values(i) + vel_values(i-1))/2;
        vel_values = sort(vel_values);
    else
        filename = strcat('opt_vel_', num2str(cputime*10000000));
        save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\velocity_results\',filename),'opt_results');
        opt_seed = x_opt_ankle; 
        i = i + 1;
    end
     damping = vel_values(i);
     
     if abs(min(diff(vel_values))) < 1e-6
         break
     end
     vel_values
end

% end
bad_stuff
too_many_iters
