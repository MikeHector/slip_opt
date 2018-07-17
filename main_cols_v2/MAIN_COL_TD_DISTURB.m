% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;

delta_TD = .002;
TDa_disturb = 0:delta_TD:.3;
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
damping = opt_results.c;

clearvars -except opt_seed apex_vel apex_height f_values delta_TD TDa_disturb damping

bad_stuff = 0;
too_many_iters = 0;
%Current TD angle
TDa_current = atan2(opt_seed(2,1), opt_seed(1,1));
TDa_values = TDa_current + TDa_disturb;

i = 1;
TDa =TDa_values(i);
bad_counter = 0;
while TDa < 1000
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(opt_seed, damping, apex_vel, apex_height, ankles_on, apex_vel, 0, TDa);
%     if (opt_results.flag <= 0) && ((damping_values(i) - damping_values(i - 1)) > 1e-3)
%         damping_values(end+1) = (damping_values(i) + damping_values(i-1))/2;
%         damping_values = sort(damping_values);           
%     else
        uniqueID = string(datetime, 'dMMyHHmmssSSSS');
        filename = strcat('opt_TDdisturb_', uniqueID);
        save(strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\TDdisturb_results\',filename),'opt_results');
%         save(strcat('D:\Documents\DRL\slip_opt\opt_results\damping_results\',filename),'opt_results');
        opt_seed = x_opt_ankle; 
        i = i + 1;
        TDa = TDa_values(i);
        
        if opt_results.flag < 0
            bad_counter = bad_counter + 1;
        elseif opt_results.flag >= 0 && bad_counter > 0
            bad_counter = 0; %reset
        end
        
        if bad_counter == 10
            break
        end
        
%     end
     disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);disp(TDa);
end

% end
bad_stuff
too_many_iters
