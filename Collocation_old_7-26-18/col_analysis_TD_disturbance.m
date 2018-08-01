% MHector
% 6.5.18
% COL analysis
clc; clear; close all
record_video = 0;
if record_video==1
    q=VideoWriter('slip_traj_over_TD_disturbance','MPEG-4');
    q.FrameRate=10;
    open(q);
end
strucc = dir('C:\Users\mike-\Google Drive\DRL- Mike Hector\SLIP Program\Collocation\optimization_results\opt_TD_disturb2*');

velmax = 3;
fig = figure;
an1 = plot(1,1); hold on
an2 = plot(2,2);
title = title('wut');
axis([-0.4, 0.3, 0, .9])
% axis([0, 0.6, -50, 50])
legend('leg torque', 'ankle torque')
xlabel('x displacement')
ylabel('y displacement')

for i = 1:length(strucc)
    filename = strucc(i).name;
    filename = strcat('C:\Users\mike-\Google Drive\DRL- Mike Hector\SLIP Program\Collocation\optimization_results\',filename);
    load(filename)
    results{i} = opt_results;
    disturb(i) = opt_results.TD_angle;
end
[disturb_sorted,i] = sort(disturb);

for k = 1:length(i)
    results_sorted_disturb{k} = results{i(k)};
end

for i = 1:length(disturb_sorted) %find(vel_sorted == vel)
    if results_sorted_disturb{i}.flag >= 1
        time = results_sorted_disturb{i}.t;
        leg_response = results_sorted_disturb{i}.Tleg;
        ankle_response = results_sorted_disturb{i}.Tankle;
        x = results_sorted_disturb{i}.x;
        y = results_sorted_disturb{i}.y;

%         an1.XData = time;
%         an1.YData = leg_response;
% 
%         an2.XData = time;
%         an2.YData = ankle_response;    
        
        an2.XData = x;
        an2.YData = y;  
        
        drawnow
        title.String = ['TD angle error = ', num2str(rad2deg(1.704 - results_sorted_disturb{i}.TD_angle))];
        pause(1)
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
