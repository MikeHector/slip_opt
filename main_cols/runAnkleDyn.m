clear; close all;
load('C:\Users\mike-\Documents\DRL\slip_opt\opt_results\opt_apex_velocity_210820182011301000.mat')
og = opt_results;
clear opt_results
modified = og;
modified.Tankle(:,1:60) = -og.Tankle;
modified.X(8,:) = modified.Tankle;
[~, new] = ankleDyn(modified);
figure; subplot(2,1,1); plot(og.x, og.y); hold on; plot(new.x, new.y)
legend('Original with ankle', 'modified ankle');
title('Max apex velocity, modified ankle torque')
xlabel('x'); ylabel('y')

subplot(2,1,2); plot(og.x, og.Tankle); hold on; plot(new.x, new.Tankle);
legend('Original with ankle', 'modified ankle');
title('Torque tapes')