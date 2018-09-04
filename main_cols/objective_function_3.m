function [ cost ] = objective_function_3( DecisionVariables, Parameters, smooth )
%objective_function This function evaluates the objective function
%   This function uses Decision variables and calculates the objective
%   function or cost
    u = DecisionVariables(7:8, :);
    T = DecisionVariables(9,1);
    v = sqrt(DecisionVariables(4,:).^2 + DecisionVariables(5,:).^2);
    r = sqrt(DecisionVariables(1,:).^2 + DecisionVariables(2,:).^2);
    hk = [DecisionVariables(9,1)/Parameters.Nstance * ones([1,Parameters.Nstance]),...
        DecisionVariables(9,2)/Parameters.Nflight * ones([1,Parameters.Nflight])];
    cost_leg = 0;
    cost_ankle = 0;
    R_leg = Parameters.R_leg;
    R_ankle = Parameters.R_ankle;
    maxXzero = MikeMax(smooth);
%     maxXzeroLeg = @(x) .2*log(1+exp(5*x));
%     maxXzeroAnkle = @(x) .05*log(1+exp(20*x));
    for i = 1:size(DecisionVariables,2)
%         if i >= Parameters.Nstance
%             disp('stop')
%         end
        cost_leg = cost_leg + ...
            R_leg * u(1, i)^2 * hk(i) +... %electrical
            maxXzero(u(1, i) * Parameters.transmission * DecisionVariables(6, i) * hk(i)); %mechanical
        cost_ankle = cost_ankle + ...
            R_ankle * u(2, i)^2 * hk(i) +... %electrical
            maxXzero(u(2, i) * Parameters.transmission_ankle *...
            (DecisionVariables(1,i) * DecisionVariables(5,i) - DecisionVariables(2,i) * DecisionVariables(4,i)) / (r(i)^2)  * hk(i)); %mechanical
    end
    cost = cost_leg + cost_ankle;
%     disp(['Cost of ankle is ', num2str(cost_ankle), ' J']);
end