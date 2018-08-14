% MHector

% 6.5.18
% COL analysis
clc; clear; close all
record_video = 0;
if record_video==1
    v=VideoWriter('Accelerating','MPEG-4');
    v.FrameRate=10;
    open(v);
end
% strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_*');%dir('D:\Documents\DRL\slip_opt\opt_results\damping_results\opt_damping_*'); %dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt_damping_*');
% strucc = dir('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\perturb_results\opt*');

% {'c', 'apex_velocity', 'disturance_f', 'TD_disturb', 'deltav', 'deltah'}
varName = 'deltav';
varmaxplot = 2.98;
plotName = 'Change in Velocity';

dirname = strcat('C:\\Users\mike-\Documents\DRL\collocation\opt_results\opt_', varName, '*');
strucc = dir(dirname);
fig = figure;
hold on
subplot(2,2,1); an1 = plot(1,1); hold on; an12 = plot(0,0,'bo');
axis([-0.3, 1, 0, 1]); title('XY Trajectory'); xlabel('X Displacement'); ylabel('Y Displacement');
subplot(2,2,2); an2 = plot(1,1); hold on; an22 = plot(2,2);
axis([-.3, 1, -13, 13]); title('Torque Trajectory'); xlabel('X Displacement'); ylabel('Torque'); legend('Ankle torque', 'Leg torque', 'Location', 'southwest')
TLmax = refline(0, 12.2); TLmax.Color = [0.8500 0.3250 0.0980]; TLmax.LineStyle = '--'; TLmax.HandleVisibility = 'off';
TLmin = refline(0, -12.2); TLmin.Color = [0.8500 0.3250 0.0980]; TLmin.LineStyle = '--'; TLmin.HandleVisibility = 'off';
TAmax = refline(0, 4.5); TAmax.Color = 'b'; TAmax.LineStyle = '--'; TAmax.HandleVisibility = 'off';
TAmin = refline(0, -4.5); TAmin.Color = 'b'; TAmin.LineStyle = '--'; TAmin.HandleVisibility = 'off';

subplot(2,2,3); an3 = plot(1,1,'ro'); hold on; an32 = plot(2,2);
axis([-.75,varmaxplot, 0, 1e3]); xlabel(plotName); ylabel('Cost');
title1 = title('wut');
subplot(2,2,4); an4 = plot(1,1);
axis([-.25, 1, -.12, .12]); xlabel('x'); ylabel('xcop')
% an2 = plot(2,2);
title('Center of Pressure')
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
        energy{q} = get_energy(results{i(k)},0);
        q = q+1;
    end
end
an32.XData = var_graph;
an32.YData = cost_graph;
% for i = 1:numel(results)
i = 1;
while results_sorted_var{i}.param.(varName) < varmaxplot
    
    if results_sorted_var{i}.param.flag > 0
        time = results_sorted_var{i}.t;
        leg_response = results_sorted_var{i}.Tleg;
        ankle_response = results_sorted_var{i}.Tankle;
        xyTraj = getXYplot(results_sorted_var{i},0);
        x = real(xyTraj.x);
        y = real(xyTraj.y);
        r = results_sorted_var{i}.r;
        k = results_sorted_var{i}.param.k;
        xcop = -ankle_response .* r ./(k .*(results_sorted_var{i}.r0 -r).* results_sorted_var{i}.y);
        
        if max(leg_response) > 12
%             pause
        end

        an1.XData = x;
        an1.YData = y;
        
        an12.XData = [results_sorted_var{i}.x(1), results_sorted_var{i}.x(end)];
        an12.YData = [results_sorted_var{i}.y(1), results_sorted_var{i}.y(end)];

        an2.XData = results_sorted_var{i}.x;
        an2.YData = ankle_response;  
        
        an22.XData = results_sorted_var{i}.x;
        an22.YData = leg_response; 
        
        an3.XData = results_sorted_var{i}.param.(varName);
        an3.YData = results_sorted_var{i}.cost;
        
        an4.XData = results_sorted_var{i}.x;
        an4.YData = xcop;

        drawnow
        title1.String = ['Energy Required when ', plotName, ' is ', num2str(results_sorted_var{i}.param.(varName)), 'm/s'];
        pause(.05)
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


legE = []; legM = []; ankleE = []; ankleM = [];
for i = 1:numel(energy)    
    legE = [legE, energy{i}.leg_e];
    legM = [legM, energy{i}.leg_m];
    ankleE = [ankleE, energy{i}.ankle_e];
    ankleM = [ankleM, energy{i}.ankle_m];
end

figure;
plot(var_graph, legE); hold on;
plot(var_graph, legM);
plot(var_graph, ankleE);
plot(var_graph, ankleM);
rl = refline(0, energy{1}.LegEtoStand); rl.LineStyle = '--';
legend('Leg electrical', 'Leg mechanical', 'Ankle electrical', 'Ankle mechanical', 'Energy to stand')