%MHector 
%7.19.18

newN = 40;
load('C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\new_obj_func\opt_damping_baseline_Rs9_last')
variableNames = {'x' 'y' 'r0' 'dx' 'dy' 'dr0' 'Tleg' 'Tankle' 't'};
oldN = opt_results.N;
oldSeed = opt_results.X;
old_opt_results = opt_results;
newSeed = zeros(length(variableNames), newN);
new_opt_results = opt_results;
new_opt_results.N = newN;

for i = 1:numel(variableNames)
    newSeed(i,:) = interp1(1:oldN, oldSeed(i,:), linspace(1,oldN,newN),'pchip');
    new_opt_results.(variableNames{i}) = newSeed(i,:);
end
new_opt_results.X = newSeed;
new_opt_results.r = sqrt(new_opt_results.x.^2 + new_opt_results.y.^2);
opt_results = new_opt_results;

