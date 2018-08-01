%MHector
%7.17.18
clc; clear;

n = 50;
m = 1; %percent

%Damping
directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt*'; %Need to have only 3 chars after \
variable_name = 'c';
%1 improvement at 55.5918

% %Velocity
% directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\velocity_results\vel*'; %Need to have only 3 chars after \
% variable_name = 'apex_velocity';

% %Delta Velocity
% directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\deltav_results\del*'; %Need to have only 3 chars after \
% variable_name = 'deltav';
% %1 improvement at 1.1263

% %Force Disturbance
% directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\fdisturb_results\opt*'; %Need to have only 3 chars after \
% variable_name = 'disturbance_f';

% %TD Disturbance
% directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\TDdisturb_results\opt*'; %Need to have only 3 chars after \
% variable_name = 'TD_disturbance';

[good, needs_improvement] = space_func(directory, variable_name, n, m);

needs_improvement


