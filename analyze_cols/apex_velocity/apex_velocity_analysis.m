%MHector 
% 8.23.18
%Let's analyze some optimization results:

%Result: at high apex velocities (2.95 m/s), the model used antagonistic
%ankle torque (pushing in the -x direction) through stance.

%Hypothesis: The slip is landing with a leg angle closer to vertical so
%that it doesn't have to spend as much energy maintaining the spring
%position. The leg motor is saturated throughout much of stance, so it's
%control authority is diminished. 

%Test: Run the optimization with the same apex velocity but with ankle
%torque off and compare the time vs height of the COM.

clear; clc; close all;
load('C:\Users\mike-\Documents\DRL\slip_opt\opt_results\opt_apex_velocity_210820182011301000.mat')
ankles = opt_results;
noAnkles = opt_results;
clear opt_results;

noAnkles.param.ankles_on = 0;
noAnkles.Tankle = zeros(size(noAnkles.Tankle));
noAnkles.X(8,:) = zeros(size(noAnkles.X(8,:)));

[~, noAnkles] = RUN_COL(noAnkles.X, noAnkles.param);

figure; 

subplot(3,1,1); plot(ankles.t, ankles.y); hold on; plot(noAnkles.t, noAnkles.y);
legend('With ankle torque', 'Without ankle torque');
title('Optimal t y trajectories'); xlabel('t'); ylabel('y');

subplot(3,1,2); plot(ankles.t, ankles.r0); hold on; plot(noAnkles.t, noAnkles.r0);
legend('With ankle torque', 'Without ankle torque');
title('r0 over t')

GRF_ankles = ankles.param.k * (ankles.r0 - ankles.r);
GRF_noAnkles = noAnkles.param.k * (noAnkles.r0 - noAnkles.r);

subplot(3,1,3); plot(ankles.t, GRF_ankles); hold on; plot(noAnkles.t, GRF_noAnkles);
legend('With ankle torque', 'Without ankle torque');
title('GRF over t')


figure; 

ankles_theta = atan2(ankles.y, ankles.x);
noAnkles_theta = atan2(noAnkles.y, noAnkles.x);

subplot(3,1,1); plot(ankles_theta, ankles.y); hold on; plot(noAnkles_theta, noAnkles.y);
legend('With ankle torque', 'Without ankle torque');
title('Optimal t y trajectories'); xlabel('t'); ylabel('y');

subplot(3,1,2); plot(ankles_theta, ankles.r0); hold on; plot(noAnkles_theta, noAnkles.r0);
legend('With ankle torque', 'Without ankle torque');
title('r0 over t')

GRF_ankles = ankles.param.k * (ankles.r0 - ankles.r);
GRF_noAnkles = noAnkles.param.k * (noAnkles.r0 - noAnkles.r);

subplot(3,1,3); plot(ankles_theta, GRF_ankles); hold on; plot(noAnkles_theta, GRF_noAnkles);
legend('With ankle torque', 'Without ankle torque');
title('GRF over t')

%Results: With torque: TD is closer to vertical... marginally. SLIP is able
%to keep it's spring at r0_start through apex. (Cost is lower with ankles).
%Stance is a little shorter with ankles. 