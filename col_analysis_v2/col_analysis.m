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
strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_*');

cmax = 410;
fig = figure;
hold on
subplot(2,2,1); an1 = plot(1,1);
axis([-0.25, .5, .25, 1]); title('xy traj'); xlabel('x'); ylabel('y');
subplot(2,2,2); an2 = plot(1,1); hold on; an22 = plot(2,2);
axis([-.25, .5, -50, 50]); title('torque traj'); xlabel('x'); ylabel('torque');
subplot(2,2,3); an3 = plot(1,1,'ro'); hold on; an32 = plot(2,2);
axis([0,cmax, 0, 10]); xlabel('damping'); ylabel('cost');
subplot(2,2,4); an4 = plot(1,1);
axis([-.25, .5, -.12, .12]); xlabel('x'); ylabel('xcop')
% an2 = plot(2,2);
title = title('wut');
% axis([-0.25, 0.25, .1, .9])

% legend('leg torque', 'ankle torque')
xlabel('time')
% ylabel('torque')

for i = 1:length(strucc)
    filename = strucc(i).name;
    filename = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\', filename);
    load(filename)
    results{i} = opt_results;
    c(i) = opt_results.c;
    if opt_results.c == 150
        pause
    end
end
[c_sorted,i] = sort(c);

q=1;
for k = 1:length(i)
    results_sorted_c{k} = results{i(k)};
    flags(k) = results{i(k)}.flag;
    if results{i(k)}.flag > 0
        c_graph(q) = results{i(k)}.c;
        cost_graph(q) = results{i(k)}.cost;
        q = q+1;
    end
end
an32.XData = c_graph;
an32.YData = cost_graph;
for i = 1:numel(results)
    if results_sorted_c{i}.flag > 0
        time = results_sorted_c{i}.t;
        leg_response = results_sorted_c{i}.Tleg;
        ankle_response = results_sorted_c{i}.Tankle;
        x = results_sorted_c{i}.x;
        y = results_sorted_c{i}.y;
        r = results_sorted_c{i}.r;
        k = results_sorted_c{i}.k;
        xcop = -ankle_response .* r ./(k .*(results_sorted_c{i}.r0 -r).* y);

        an1.XData = x;
        an1.YData = y;

        an2.XData = x;
        an2.YData = ankle_response;  
        
        an22.XData = x;
        an22.YData = leg_response; 
        
        an3.XData = results_sorted_c{i}.c;
        an3.YData = results_sorted_c{i}.cost;
        
        an4.XData = x;
        an4.YData = xcop;

        drawnow
        title.String = ['damping = ', num2str(floor(results_sorted_c{i}.c))];
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
