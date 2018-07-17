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
strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\TDdisturb_results\opt*'); %DRL desktop

disturbMax = 150;
fig = figure;
hold on
subplot(2,2,1); an1 = plot(1,1);
axis([-0.3, .3, .25, 1]); title('xy traj'); xlabel('x'); ylabel('y');
subplot(2,2,2); an2 = plot(1,1); hold on; an22 = plot(2,2);
axis([-.3, .3, -5, 50]); title('torque traj'); xlabel('x'); ylabel('torque');
subplot(2,2,3); an3 = plot(1,1,'ro'); hold on; an32 = plot(2,2);
axis([0,2, 0, .5]); xlabel('TD Disturbance'); ylabel('cost');
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
    filename = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\TDdisturb_results\', filename); %DRL Desktop
    load(filename)
    results{i} = opt_results;
    TDa(i) = atan2(opt_results.y(1),opt_results.x(1));
    if opt_results.disturbance_f == 50
%         pause
    end
end
[TDa_sorted,i] = sort(TDa);

q=1;
for k = 1:length(i)
    results_sorted_TD{k} = results{i(k)};
    flags(k) = results{i(k)}.flag;
    if results{i(k)}.flag > 0
        TD_graph(q) = results{i(k)}.disturbance_f;
        cost_graph(q) = results{i(k)}.cost;
        q = q+1;
    end
end
an32.XData = TD_graph;
an32.YData = cost_graph;

for i = 1:numel(results)
    %Make energy shared figure
    energy_leg(i) = sum(results_sorted_TD{i}.Tleg.^2);
    energy_ankle(i) = sum(results_sorted_TD{i}.Tankle.^2);
    
    if results_sorted_TD{i}.flag > 0
        time = results_sorted_TD{i}.t;
        leg_response = results_sorted_TD{i}.Tleg;
        ankle_response = results_sorted_TD{i}.Tankle;
        x = results_sorted_TD{i}.x;
        y = results_sorted_TD{i}.y;
        r = results_sorted_TD{i}.r;
        k = results_sorted_TD{i}.k;
        xcop = -ankle_response .* r ./(k .*(results_sorted_TD{i}.r0 -r).* y);
        
        TDangle = atan2(results_sorted_TD{i}.y(1), results_sorted_TD{i}.x(1));
        
        if max(leg_response) > 49.9
%             pause
        end

        an1.XData = x;
        an1.YData = y;

        an2.XData = x;
        an2.YData = ankle_response;  
        
        an22.XData = x;
        an22.YData = leg_response; 
        
        an3.XData = TDangle;
        an3.YData = results_sorted_TD{i}.cost;
        
        an4.XData = x;
        an4.YData = xcop;

        drawnow
        title1.String = ['TD angle= ', num2str(TDangle)];
        pause(.5)
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
[TDUnique, indUnique] = unique(TDa);
% barArray = [energy_leg(indUnique)' energy_ankle(indUnique)']; 
% figure
% abar = bar(cUnique',barArray, 'stacked');
figure;
plot(TDUnique, energy_leg(indUnique)); hold on; 
plot(TDUnique, energy_ankle(indUnique),'r');
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
