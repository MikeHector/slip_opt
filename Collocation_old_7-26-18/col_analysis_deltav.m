% MHector
% 6.5.18
% COL analysis
clc; clear; close all
record_video = 0;
if record_video==1
    q=VideoWriter('torques_over_deltav','MPEG-4');
    q.FrameRate=10;
    open(q);
end
strucc = dir('opt_deltav_*');

velmax = 3;
fig = figure;
an1 = plot(1,1); hold on
an2 = plot(2,2);
title = title('wut');
% axis([-0.25, 0.25, .1, .9])
legend('leg torque', 'ankle torque')
xlabel('time')
ylabel('torque')

for i = 1:length(strucc)
    filename = strucc(i).name;
    load(filename)
    results{i} = opt_results;
    deltav(i) = opt_results.end_vel - opt_results.apex_velocity;
end
[deltav_sorted,i] = sort(deltav);

for k = 1:length(i)
    results_sorted_vel{k} = results{i(k)};
end

for i = 1:length(deltav_sorted) %find(vel_sorted == vel)
    [results_sorted_vel{i}.end_vel - results_sorted_vel{i}.apex_velocity, results_sorted_vel{i}.flag]
    pause(.5)
    if results_sorted_vel{i}.flag >= 0
        time = results_sorted_vel{i}.t;
        leg_response = results_sorted_vel{i}.Tleg;
        ankle_response = results_sorted_vel{i}.Tankle;
        x = results_sorted_vel{i}.x;
        y = results_sorted_vel{i}.y;

        an1.XData = time;
        an1.YData = leg_response;

        an2.XData = time;
        an2.YData = ankle_response;    
        
%         an2.XData = x;
%         an2.YData = y;  
        
        drawnow
        title.String = ['deltav = ', num2str(results_sorted_vel{i}.end_vel - results_sorted_vel{i}.apex_velocity)];
        pause(.05)
        if record_video==1
            F=getframe(gcf);
            writeVideo(q,F);
        end
    end
end

if record_video == 1
    close(q)
end

% figure
% subplot(2,2,1)
% 
% plot(c,flag,'bo')
% title('fmincon ending state flag')
% subplot(2,2,2)
% 
% plot(c,cost,'bo')
% title('cost')
% subplot(2,2,3)
% 
% plot(c,leg_response,'bo'); hold on
% plot(c,ankle_response,'ro')
% title('torque squared')
