function [ cost ] = objective_function( DecisionVariables, Parameters )
%objective_function This function evaluates the objective function
%   This function uses Decision variables and calculates the objective
%   function or cost
    u = DecisionVariables(7:8, :);
    T = DecisionVariables(9,1);
    hk = T / Parameters.N;
    cost_leg = 0;
    cost_ankle = 0;
    for i = 1:Parameters.N - 1
        cost_leg = cost_leg + (hk / 2 * (u(1, i)^2 + u(1, i + 1)^2)) * hk;
        cost_ankle = cost_ankle + (hk / 2 * (u(2, i)^2 + u(2, i + 1)^2)) * hk;
    end
    cost = cost_leg + cost_ankle;
end