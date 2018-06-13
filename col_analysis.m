% MHector
% 6.5.18
% COL analysis
clc; clear; close all
record_video = 0;
if record_video==1
    v=VideoWriter('torques_over_damping','MPEG-4');
    v.FrameRate=10;
    open(v);
end
strucc = dir('opt_damping_*');

cmax = 450;
fig = figure;
an1 = plot(1,1); hold on
% an2 = plot(2,2);
title = title('wut');
axis([-0.25, 0.25, .1, .9])
% legend('leg torque', 'ankle torque')
xlabel('time')
% ylabel('torque')

for i = 1:length(strucc)
    filename = strucc(i).name;
    load(filename)
    results{i} = opt_results;
    c(i) = opt_results.c;
end
[c_sorted,i] = sort(c);

for k = 1:length(i)
    results_sorted_c{k} = results{i(k)};
end

for i = 1:find(c_sorted == cmax)
    if results_sorted_c{i}.flag > 0
        time = results_sorted_c{i}.t;
        leg_response = results_sorted_c{i}.Tleg;
        ankle_response = results_sorted_c{i}.Tankle;
        x = results_sorted_c{i}.x;
        y = results_sorted_c{i}.y;

        an1.XData = x;
        an1.YData = y;

%         an2.XData = time;
%         an2.YData = ankle_response;    
        
        drawnow
        title.String = ['damping = ', num2str(results_sorted_c{i}.c)];
        pause(.05)
        if record_video==1
            F=getframe(gcf);
            writeVideo(v,F);
        end
    end
end

if record_video == 1
    close(v)
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
