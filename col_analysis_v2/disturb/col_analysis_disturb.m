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
% strucc = dir('D:\Documents\DRL\slip_opt\opt_results\velocity_results\vel*');  %My desktop
strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\fdisturb_results\opt*'); %DRL desktop

disturbMax = 50;
fig = figure;
hold on
subplot(2,2,1); an1 = plot(1,1);
axis([-0.3, .3, .25, 1]); title('xy traj'); xlabel('x'); ylabel('y');
subplot(2,2,2); an2 = plot(1,1); hold on; an22 = plot(2,2);
axis([-.3, .3, -5, 50]); title('torque traj'); xlabel('x'); ylabel('torque');
subplot(2,2,3); an3 = plot(1,1,'ro'); hold on; an32 = plot(2,2);
axis([0,disturbMax, 0, .5]); xlabel('velocity'); ylabel('cost');
subplot(2,2,4); an4 = plot(1,1);
axis([-.25, .5, -.12, .12]); xlabel('x'); ylabel('xcop')
% an2 = plot(2,2);
title1 = title('wut');
% axis([-0.25, 0.25, .1, .9])

% legend('leg torque', 'ankle torque')
xlabel('time')
% ylabel('torque')

for i = 1:length(strucc)
    filename = strucc(i).name;
%     filename = strcat('D:\Documents\DRL\slip_opt\opt_results\velocity_results\', filename); %My desktop
    filename = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\fdisturb_results\', filename); %DRL Desktop
    load(filename)
    results{i} = opt_results;
    disturb(i) = opt_results.disturbance_f;
    if opt_results.disturbance_f == 50
        pause
    end
end
[disturb_sorted,i] = sort(disturb);

q=1;
for k = 1:length(i)
    results_sorted_disturb{k} = results{i(k)};
    flags(k) = results{i(k)}.flag;
    if results{i(k)}.flag > 0
        disturb_graph(q) = results{i(k)}.disturbance_f;
        cost_graph(q) = results{i(k)}.cost;
        q = q+1;
    end
end
an32.XData = disturb_graph;
an32.YData = cost_graph;

for i = 1:numel(results)
    %Make energy shared figure
    energy_leg(i) = sum(results_sorted_disturb{i}.Tleg.^2);
    energy_ankle(i) = sum(results_sorted_disturb{i}.Tankle.^2);
    
    if results_sorted_disturb{i}.flag > 0
        time = results_sorted_disturb{i}.t;
        leg_response = results_sorted_disturb{i}.Tleg;
        ankle_response = results_sorted_disturb{i}.Tankle;
        x = results_sorted_disturb{i}.x;
        y = results_sorted_disturb{i}.y;
        r = results_sorted_disturb{i}.r;
        k = results_sorted_disturb{i}.k;
        xcop = -ankle_response .* r ./(k .*(results_sorted_disturb{i}.r0 -r).* y);
        
        if max(leg_response) > 49.9
%             pause
        end

        an1.XData = x;
        an1.YData = y;

        an2.XData = x;
        an2.YData = ankle_response;  
        
        an22.XData = x;
        an22.YData = leg_response; 
        
        an3.XData = results_sorted_disturb{i}.disturbance_f;
        an3.YData = results_sorted_disturb{i}.cost;
        
        an4.XData = x;
        an4.YData = xcop;

        drawnow
        title1.String = ['disturbance force = ', num2str(results_sorted_disturb{i}.disturbance_f)];
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

%Make energy shared figure
[disturbUnique, indUnique] = unique(disturb);
% barArray = [energy_leg(indUnique)' energy_ankle(indUnique)']; 
% figure
% abar = bar(cUnique',barArray, 'stacked');
figure;
plot(disturbUnique, energy_leg(indUnique)); hold on; 
plot(disturbUnique, energy_ankle(indUnique),'r');
% a = line([1 1],[0, 12*10^4]); a.LineStyle = '--';
xlabel('Apex Velocity')
ylabel('Energy')
legend('Leg energy', 'Ankle energy')
title('Optimal energies of actuators through stance')
% a.Color = 'k';
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
