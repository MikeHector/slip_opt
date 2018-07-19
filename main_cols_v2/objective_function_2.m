function [ cost ] = objective_function_2( DecisionVariables, Parameters )
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
    R = Parameters.R;
    for i = 1:Parameters.N - 1
        cost_leg = cost_leg + ...
            R * u(1, i)^2 * hk + max(0, u(1, i) * Parameters.transmission * DecisionVariables(6, i) * hk);
        cost_ankle = cost_ankle + ...
            R * u(2, i)^2 * hk + max(0, u(2, i) / r(i) * v(i) * hk);
    end
    cost = cost_leg + cost_ankle;
end