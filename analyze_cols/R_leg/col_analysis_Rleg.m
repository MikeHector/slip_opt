% MHector
% 8.14.18
% COL analysis
clc; clear; close all
record_video = 0;
if record_video==1
    v=VideoWriter('Damping','MPEG-4');
    v.FrameRate=10;
    open(v);
end

% {'c', 'apex_velocity', 'disturance_f', 'TD_disturb', 'deltav', 'deltah'}
varName = 'R_leg';
varmaxplot = 25;
varminplot = 0;
energyMax = 800;
plotName = 'R leg';
cf = pwd; %Path stuff
addpath(strcat(cf(1:strfind(pwd, 'slip_opt')-1), 'slip_opt\main_cols\')); %Add main col folder to path
dirComp = getSaveDir('Michael-PC'); %Change if you're running on a different computer

dirname = strcat(dirComp, 'opt_', varName, '*');
strucc = dir(dirname);
fig = figure;
hold on
subplot(2,2,1); an1 = plot(1,1);
axis([0,1, 0, 1.2]); title('Y Height Through Cycle'); xlabel('Normalized Time'); ylabel('Y Displacement');
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
axis([0,1, 0, 4]); xlabel('Normalized Time'); ylabel('GRF [N]')
% an2 = plot(2,2);
title('GRF normalized by weight')
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
    if varr(i) > 1.1 && varr(i) < 1.3
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
        energy{q} = get_energy2(results{i(k)},0);
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
        an4.YData = results_sorted_var{i}.param.k * (results_sorted_var{i}.r0 - results_sorted_var{i}.r)/...
                    (results_sorted_var{i}.param.m * results_sorted_var{i}.param.g);

        drawnow
        title1.String = ['Energy Required when ', plotName, ' is ', num2str(results_sorted_var{i}.param.(varName)), 'm/s'];
%         fease = results_sorted_var{i}.param.fmincon_stuff.constrviolation;
%         disp(['constraint violation is ', num2str(fease)])
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
        title2.String = ['Ground Reaction Forces when ', plotName, ' is ', num2str(results_sorted_var{i}.param.(varName))];
        
        leglen.XData = time;
        leglen.YData = r;
        drawnow;
        pause(.005);
   end
    i = i + 1;
end
opacityscale = linspace(0,1,i);
opacityscaleflip = flip(opacityscale);
figure;
i = 1;
while results_sorted_var{i}.param.(varName) < varmaxplot
    
   if results_sorted_var{i}.param.flag > 0
        time = results_sorted_var{i}.t/results_sorted_var{i}.t(end);
        r = results_sorted_var{i}.r;
        k = results_sorted_var{i}.param.k;
        GRF = k* (results_sorted_var{i}.r0 - results_sorted_var{i}.r)/ (results_sorted_var{i}.param.m * results_sorted_var{i}.param.g);
        
    subplot(3,1,1)
    hold on;
    plot(time, GRF,'Color', [opacityscale(i),opacityscaleflip(i),opacityscaleflip(i)])
%     alpha(opacityscale)
    subplot(3,1,2)
    hold on;
    plot(time, r,'Color', [opacityscale(i),opacityscaleflip(i),opacityscaleflip(i)])
    subplot(3,1,3)
    hold on;
    plot(var_graph(i), cost_graph(i), 'o', 'Color', [opacityscale(i),opacityscaleflip(i),opacityscaleflip(i)])
%     alpha(opacityscale)
%         title2.String = ['Ground Reaction Forces when ', plotName, ' is ', num2str(results_sorted_var{i}.param.(varName))];
        
%         leglen.XData = time;
%         leglen.YData = r;
        drawnow;
        pause(.05);
   end
    i = i + 1;
end
subplot(3,1,1); grf = plot(0,0); xlabel('Normalized Time'); ylabel('GRF Normalized by Weight'); title('Ground Reaction Force');
subplot(3,1,2); leglen = plot(0,0); xlabel('Normalized Time'); ylabel('Leg Length'); title('Leg Length Through Stance')
subplot(3,1,3); xlabel('Electrical losses factor'); ylabel('Energy required'); title('Cost for equilibrium gait cycle')
hold on; plot(1.4, 32.69, 'r*'); plot(1.4, 32.69, 'ro'); text(.5, 100, 'Cassie losses factor')