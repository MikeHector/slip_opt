function [ cost ] = objective_function_3( DecisionVariables, Parameters )
%objective_function This function evaluates the objective function
%   This function uses Decision variables and calculates the objective
%   function or cost
    u = DecisionVariables(7:8, :);
    T = DecisionVariables(9,1);
    v = sqrt(DecisionVariables(4,:).^2 + DecisionVariables(5,:).^2);
    r = sqrt(DecisionVariables(1,:).^2 + DecisionVariables(2,:).^2);
    hk = T / Parameters.N;
    cost_leg = 0;
    cost_ankle = 0;
    R_leg = Parameters.R_leg;
    R_ankle = Parameters.R_ankle;
    for i = 1:Parameters.N
        cost_leg = cost_leg + ...
            R_leg * u(1, i)^2 * hk +... %electrical
            max(0, u(1, i) * Parameters.transmission * DecisionVariables(6, i) * hk); %mechanical
        cost_ankle = cost_ankle + ...
            R_ankle * u(2, i)^2 * hk +... %electrical
            max(0, u(2, i) * Parameters.transmission_ankle *...
            (DecisionVariables(1,i) * DecisionVariables(5,i) - DecisionVariables(2,i) * DecisionVariables(4,i)) / (r(i)^2)  * hk); %mechanical
    end
    cost = cost_leg + cost_ankle;
end