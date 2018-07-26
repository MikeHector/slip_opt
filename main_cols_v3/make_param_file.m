%Import old opt_results
%for varName = {'c', 'deltav', 'FDisturb', 'TDA', 'apexVel'}
clear
load('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\new_obj_func\opt_damping_baseline_Rs9_short_newerDyn_faster_higher.mat')
fieldNames = fieldnames(opt_results);
for varInd = 1:numel(fieldNames) - 14
    var = fieldNames{varInd};
    param.(var) = opt_results.(var);
end
bad_fields = {'end_vel', 'apex_height_final', 'obj_func', 'y_start', 'x_start',...
    'x_dot_start', 'y_dot_start'};
rmfield(param, bad_fields);
param.deltav = 0;
param.deltah = 0; 
rmfield(param, bad_fields) %Why do I need to do this twice?