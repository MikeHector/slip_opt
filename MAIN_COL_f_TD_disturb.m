% MHector
% 6.1.18
% Master script to run a bunch of collocations from

%Let's run over several damping values using previous damping value
%solution as a seed for the next optimization
clear; clc;


% %Load initial seed
% load('C:\Users\mike-\Google Drive\DRL- Mike Hector\SLIP Program\Collocation\SLIP Collocation_v3\opt_f_disturb_baseline') %Baseline from other analysis
% seed = opt_results.X;
% apex_vel = opt_results.apex_velocity;
% end_vel = opt_results.apex_velocity;
% apex_height = opt_results.apex_height;
% damping_values = opt_results.c;
% disturbance_f = 0;
% baseline_TD = atan2(opt_results.y(1), opt_results.x(1));
% TD_disturb = linspace(0,.3,50);
% TD_angle = baseline_TD + TD_disturb;
% clear opt_results
% %Save to my googledrive
% savelocation = 'C:\Users\mike-\Google Drive\DRL- Mike Hector\SLIP Program\Collocation\optimization_results\';
% bad_stuff = 0;
% too_many_iters = 0;
% % wb = waitbar(0, 'Optimizing unoptimized optimizations');
% %Start from seed and go up
% for i = 1:length(TD_angle)
%     ankles_on = 1;
%     [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel, apex_height, ankles_on, end_vel, disturbance_f, TD_angle(i));
%     fileLocName = strcat(savelocation,'opt_TD_disturb2_',num2str(cputime*10000000));
%     save(fileLocName,'opt_results');
%     seed = x_opt_ankle;
%     if opt_results.flag < 0
%         bad_stuff = bad_stuff + 1;
%     elseif opt_results.flag == 0
%         too_many_iters = too_many_iters + 1;
%     end
% %     close(wb)
% %     wb = waitbar(i/length(TD_angle),'Optimizing unoptimized optimizations');
% end
% too_many_iters
% bad_stuff
% clear

too_many_iters = 0;
bad_stuff = 0 ;
%Load initial seed
load('C:\Users\mike-\Google Drive\DRL- Mike Hector\SLIP Program\Collocation\SLIP Collocation_v3\opt_f_disturb_baseline') %Baseline from other analysis
seed = opt_results.X;
apex_vel = opt_results.apex_velocity;
end_vel = opt_results.apex_velocity;
apex_height = opt_results.apex_height;
damping_values = opt_results.c;
disturbance_f = 0;
baseline_TD = atan2(opt_results.y(1), opt_results.x(1));
TD_disturb = linspace(0,-.3,50);
TD_angle = baseline_TD + TD_disturb;
clear opt_results
%Save to my googledrive
savelocation = 'C:\Users\mike-\Google Drive\DRL- Mike Hector\SLIP Program\Collocation\optimization_results\';


% wb = waitbar(0, 'Optimizing unoptimized optimizations');
%Start from seed and go up
for i = 1:length(TD_angle)
    ankles_on = 1;
    [x_opt_ankle, opt_results] = RUN_COL(seed, damping_values, apex_vel, apex_height, ankles_on, end_vel, disturbance_f, TD_angle(i));
    fileLocName = strcat(savelocation,'opt_TD_disturb2_',num2str(cputime*10000000));
    save(fileLocName,'opt_results');
    seed = x_opt_ankle;
    if opt_results.flag < 0
        bad_stuff = bad_stuff + 1;
    elseif opt_results.flag == 0
        too_many_iters = too_many_iters + 1;
    end
%     close(wb)
%     wb = waitbar(i/length(TD_angle),'Optimizing unoptimized optimizations');
end

bad_stuff
too_many_iters
