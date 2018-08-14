% MHector

% 6.5.18
% COL analysis
clc; clear; close all
record_video = 0;
if record_video==1
    v=VideoWriter('Final Apex Heights','MPEG-4');
    v.FrameRate=10;
    open(v);
end
% strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_*');%dir('D:\Documents\DRL\slip_opt\opt_results\damping_results\opt_damping_*'); %dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_*');
% strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\perturb_results\opt*');

% {'c', 'apex_velocity', 'disturance_f', 'TD_disturb', 'deltav', 'deltah'}
varName = 'deltah';
varmaxplot = .4899;
varminplot = -.1;
plotName = 'Final Apex Height';

dirname = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\opt_', varName, '*');
strucc = dir(dirname);
assert(numel(strucc) > 0, 'No files by that name in the directory')
fig = figure;
hold on

%XY Trajectory
subplot(2,2,1); an1 = plot(1,1);
axis([-0.3, .3, 0, 1.5]); title('XY Trajectory'); xlabel('X Displacement'); ylabel('Y Displacement');

%Torques
subplot(2,2,2); an2 = plot(1,1); hold on; an22 = plot(2,2);
axis([-.3, .3, -13, 13]); title('Torque Trajectory'); xlabel('X Displacement'); ylabel('Torque'); legend('Ankle torque', 'Leg torque', 'Location', 'southwest')
TLmax = refline(0, 12.2); TLmax.Color = [0.8500 0.3250 0.0980]; TLmax.LineStyle = '--'; TLmax.HandleVisibility = 'off';
TLmin = refline(0, -12.2); TLmin.Color = [0.8500 0.3250 0.0980]; TLmin.LineStyle = '--'; TLmin.HandleVisibility = 'off';
TAmax = refline(0, 4.5); TAmax.Color = 'b'; TAmax.LineStyle = '--'; TAmax.HandleVisibility = 'off';
TAmin = refline(0, -4.5); TAmin.Color = 'b'; TAmin.LineStyle = '--'; TAmin.HandleVisibility = 'off';

subplot(2,2,3); an3 = plot(1,1,'ro'); hold on; an32 = plot(2,2);
axis([varminplot,varmaxplot, 0, 850]); xlabel('Final Apex Height'); ylabel('Cost');
title1 = title('wut');

subplot(2,2,4); an4 = plot(1,1);
axis([-.3, .3, -.12, .12]); xlabel('X Displacement'); ylabel('Center of Pressure')
title('Center of Pressure');

% an2 = plot(2,2);

% axis([-0.25, 0.25, .1, .9])

% legend('leg torque', 'ankle torque')
% xlabel('time')
% ylabel('torque')

for i = 1:length(strucc)
    filename = strucc(i).name;
    filename = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\', filename); %strcat('D:\Documents\DRL\slip_opt\opt_results\damping_results\', filename); 
    load(filename)
    results{i} = opt_results;
    varr(i) = opt_results.param.(varName);
    if varr(i) == 0
%         pause
    end
end
[var_sorted,i] = sort(varr);

q=1;
for k = 1:length(i)
    results_sorted_var{k} = results{i(k)};
    flags(k) = results{i(k)}.param.flag;
    if results{i(k)}.param.flag > 0
        var_graph(q) = results{i(k)}.param.(varName);
        cost_graph(q) = results{i(k)}.cost;
        q = q+1;
    end
    max(results_sorted_var{k}.Tankle)
end
an32.XData = var_graph;
an32.YData = cost_graph;
% for i = 1:numel(results)
i = 1;
while results_sorted_var{i}.param.(varName) < varmaxplot
    %Make energy shared figure
    energy_leg(i) = sum(results_sorted_var{i}.Tleg.^2);
    energy_ankle(i) = sum(results_sorted_var{i}.Tankle.^2);
    
    if results_sorted_var{i}.param.flag > 0 && results_sorted_var{i}.param.(varName) > varminplot
        time = results_sorted_var{i}.t;
        leg_response = results_sorted_var{i}.Tleg;
        ankle_response = results_sorted_var{i}.Tankle;
        x = results_sorted_var{i}.x;
        y = results_sorted_var{i}.y;
        r = results_sorted_var{i}.r;
        k = results_sorted_var{i}.param.k;
        xcop = -ankle_response .* r ./(k .*(results_sorted_var{i}.r0 -r).* y);
        
        if max(leg_response) > 12
%             pause
        end

        an1.XData = x;
        an1.YData = y;

        an2.XData = x;
        an2.YData = ankle_response;  
        
        an22.XData = x;
        an22.YData = leg_response; 
        
        an3.XData = results_sorted_var{i}.param.(varName);
        an3.YData = results_sorted_var{i}.cost;
        
        an4.XData = x;
        an4.YData = xcop;

        drawnow
        title1.String = ['Energy Required when ', plotName, ' is ', num2str(1+results_sorted_var{i}.param.(varName)), 'm'];
        pause(.1)
        if record_video==1
            F=getframe(gcf);
            writeVideo(v,F);
        end
    end
    i = i + 1;
end

if record_video == 1
    close(v)
end

%Make energy shared figure
% [varUnique, indUnique] = unique(var);
% % barArray = [energy_leg(indUnique)' energy_ankle(indUnique)']; 
% % figure
% % abar = bar(cUnique',barArray, 'stacked');
% figure;
% plot(varUnique, energy_leg(indUnique)); hold on; 
% plot(varUnique, energy_ankle(indUnique),'r');
% % a = line([71 71],[0, 12*10^4]); a.LineStyle = '--';
% xlabel(plotName)
% ylabel('Energy')
% legend('Leg energy', 'Ankle energy','Leg begins saturation')
% title('Optimal energies of actuators through stance')
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
