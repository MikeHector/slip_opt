function [ cost ] = objective_function_3( DecisionVariables, Parameters )
%objective_function This function evaluates the objective function
%   This function uses Decision variables and calculates the objective
%   function or cost
    x = DecisionVariables(1, :);
%     x(43);
    y = DecisionVariables(2, :);
    r0 = DecisionVariables(3, :);
    dx = DecisionVariables(4, :);
    dy = DecisionVariables(5, :);
    dr0 = DecisionVariables(6, :);
    Tleg = DecisionVariables(7, :);
    Tankle = DecisionVariables(8, :);
    r = sqrt(x.^2 + y.^2);
    Tstance = DecisionVariables(9,1);
    Tflight = DecisionVariables(9,2);
    hk = [Tstance/Parameters.Nstance * ones([1,Parameters.Nstance]),...
        Tflight/Parameters.Nflight * ones([1,Parameters.Nflight])];
    cost_leg = 0;
    cost_ankle = 0;
    R_leg = Parameters.R_leg;
    R_ankle = Parameters.R_ankle;
    maxXzero = MikeMax(1);
    for i = 1:size(DecisionVariables,2)
        cost_leg = cost_leg + ...
            R_leg * Tleg(i)^2 * hk(i) +... %electrical
            maxXzero(Tleg(i) * Parameters.transmission * dr0(i) * hk(i)); %mechanical
        cost_ankle = cost_ankle + ...
            R_ankle * Tankle(i)^2 * hk(i) +... %electrical
            maxXzero(Tankle(i) * Parameters.transmission_ankle *...
            (x(i) * dy(i) - y(i) * dx(i)) / (r(i)^2)  * hk(i)); %mechanical
    end
    cost = cost_leg + cost_ankle;
end