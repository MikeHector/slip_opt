%Import old opt_results
%for varName = {'c', 'deltav', 'FDisturb', 'TDA', 'apexVel'}
clear
load('baseline')
% fieldNames = fieldnames(opt_results.param);
goodfields = {'m', 'k', 'r0_start', 'g', 'anklemax', 'legmax', 'N', 'dof', 'cntrl_dof'...
    'i_motor', 'transmission', 'lf', 'r0_min', 'r0_max', 'timemax', 'timemin', 'c',...
    'fmincon_stuff', 'disturbance_f', 'disturbance_t_start',...
    'disturbance_t_end', 'R_leg', 'R_ankle', 'mechA_leg', 'mechA_ankle',...
    'transmission_ankle', 'apex_height', 'apex_velocity', 'options', 'deltav',...
    'deltah', 'ankles_on'};
for varInd = 1:numel(goodfields)
    param.(goodfields{varInd}) = opt_results.param.(goodfields{varInd});
end

param.lf = .165;
param.baseline_TDA = atan2(opt_results.y(1),opt_results.x(1));
param.TD_disturb = NaN;
R.param = param;

keepField = {'x', 'y', 'r0', 'dx', 'dy', 'dr0', 'Tleg', 'Tankle', 'T', 't', 'r', 'cost', 'X'};
for fieldInd = 1:numel(keepField)
    R.(keepField{fieldInd}) = opt_results.(keepField{fieldInd});
end

