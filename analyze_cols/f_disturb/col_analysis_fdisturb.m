% MHector
% 8.14.18
% COL analysis
clc; clear; close all
record_video = 1;
if record_video==1
    v=VideoWriter('Accelerating','MPEG-4');
    v.FrameRate=10;
    open(v);
end

% {'c', 'apex_velocity', 'disturbance_f', 'TD_disturb', 'deltav', 'deltah'}
varName = 'disturbance_f';
varmaxplot = 100;
varminplot = -100;
energyMax = 1000;
plotName = 'Force Disturbance During Stance';
cf = pwd; %Path stuff
addpath(strcat(cf(1:strfind(pwd, 'collocation')-1), 'collocation\main_cols\')); %Add main col folder to path
dirComp = getSaveDir('DRL-PC'); %Change if you're running on a different computer

dirname = strcat(dirComp, 'opt_', varName, '*');
strucc = dir(dirname);
fig = figure;
hold on
subplot(2,2,1); an1 = plot(1,1);
axis([0,1, 0, 1.5]); title('Y Height Through Cycle'); xlabel('Normalized Time'); ylabel('Y Displacement');
subplot(2,2,2); an2 = plot(1,1); hold on; an22 = plot(2,2);
axis([0, 1, -13, 13]); title('Torque Trajectory'); xlabel('Normalized Time'); ylabel('Torque'); legend('Ankle torque', 'Leg torque', 'Location', 'southwest')
TLmax = refline(0, 12.2); TLmax.Color = [0.8500 0.3250 0.0980]; TLmax.LineStyle = '--'; TLmax.HandleVisibility = 'off';
TLmin = refline(0, -12.2); TLmin.Color = [0.8500 0.3250 0.0980]; TLmin.LineStyle = '--'; TLmin.HandleVisibility = 'off';
TAmax = refline(0, 4.5); TAmax.Color = 'b'; TAmax.LineStyle = '--'; TAmax.HandleVisibility = 'off';
TAmin = refline(0, -4.5); TAmin.Color = 'b'; TAmin.LineStyle = '--'; TAmin.HandleVisibility = 'off';

subplot(2,2,3); an3 = plot(1,1,'ro'); hold on; an32 = plot(2,2);
axis([varminplot,varmaxplot, 0, energyMax]); xlabel(plotName); ylabel('Cost');
title1 = title('wut');
subplot(2,2,4); an4 = plot(1,1);
axis([0,1, -.12, .12]); xlabel('Normalized Time'); ylabel('xcop')
% an2 = plot(2,2);
title('Center of Pressure')
% axis([-0.25, 0.25, .1, .9])

% legend('leg torque', 'ankle torque')
% xlabel('time')
% ylabel('torque')

for i = 1:length(strucc)
    filename = strucc(i).name;
    filename = strcat(dirComp, filename);
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
        time = results_sorted_var{i}.t / results_sorted_var{i}.t(end);
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

        an1.XData = real(xyTraj.t)/real(xyTraj.t(end));
        an1.YData = y;
        
        an12.XData = time;
        an12.YData = [results_sorted_var{i}.y(1), results_sorted_var{i}.y(end)];

        an2.XData = time;
        an2.YData = ankle_response;  
        
        an22.XData = time;
        an22.YData = leg_response; 
        
        an3.XData = results_sorted_var{i}.param.(varName);
        an3.YData = results_sorted_var{i}.cost;
        
        an4.XData = time;
        an4.YData = xcop;

        drawnow
        title1.String = ['Energy Required when ', plotName, ' is ', num2str(results_sorted_var{i}.param.(varName)), ' N'];
        pause(.005)
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

% Energy figure

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

% GRF/leg length graph

figure;
subplot(2,1,1); grf = plot(0,0); xlabel('Normalized Time'); ylabel('GRF Normalized by Weight'); title2 = title('Ground Reaction Force');
axis([0 1 0 3.2])
subplot(2,1,2); leglen = plot(0,0); xlabel('Normalized Time'); ylabel('Leg Length'); title('Leg Length Through Stance')
axis([0 1 0 1.2])

i = 1;
while results_sorted_var{i}.param.(varName) < varmaxplot
    
   if results_sorted_var{i}.param.flag > 0
        time = results_sorted_var{i}.t/results_sorted_var{i}.t(end);
        r = results_sorted_var{i}.r;
        k = results_sorted_var{i}.param.k;
        GRF = k* (results_sorted_var{i}.r0 - results_sorted_var{i}.r)/ (results_sorted_var{i}.param.m * results_sorted_var{i}.param.g);
        
        grf.XData = time;
        grf.YData = GRF;
        title2.String = ['Ground Reaction Forces when ', plotName, ' is ', num2str(results_sorted_var{i}.param.(varName)),' N'];
        
        leglen.XData = time;
        leglen.YData = r;
        drawnow;
        pause(.05);
   end
    i = i + 1;
end

