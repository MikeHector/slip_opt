% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

delta_v = .01;
end_vel_values = 0:delta_damping:2;
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

load('C://Users/mike-/Documents/DRL/collocation/main_cols_v2/opt_damping_150.mat')
opt_seed = opt_results.X;
damping = opt_results.c;
clear opt_results

bad_stuff = 0;
too_many_iters = 0;

i = 1;
badcounter = 0;
for i = 1:length(end_vel_values)
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel +end_vel_values(i), 0, atan2(opt_seed(2,1),opt_seed(1,1)));
%     if opt_results.flag <= 0
%         damping_values(end+1) = (damping_values(i) + damping_values(i-1))/2;
%         damping_values = sort(damping_values);
%         if abs(damping_values(i) - damping_values(i - 1)) < 1e-6
%             i = i + 1;
%         end
%     else
        filename = strcat('opt_damping_', num2str(cputime*10000000));
        save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\',filename),'opt_results');
        opt_seed = x_opt_ankle; 
%         i = i + 1;
%     end
    
    if opt_results.flag < 0
        badcounter = badcounter + 1;
    elseif opt_results.flag > 0 && badcounter > 0
        badcounter = 0;
    end    
     
     if opt_results.flag == 0
         too_many_iters = too_many_iters + 1;
     elseif opt_results.flag < 0
         bad_stuff = bad_stuff + 1;
     end
     
     if badcounter == 5
         break
     end
end

% end
bad_stuff
too_many_iters
