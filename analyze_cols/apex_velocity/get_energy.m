function [cost] = get_energy(OPT_RES, plotifone)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    u = OPT_RES.X(7:8, :);
    T = OPT_RES.X(9,1);
    v = sqrt(OPT_RES.X(4,:).^2 + OPT_RES.X(5,:).^2);
    r = sqrt(OPT_RES.X(1,:).^2 + OPT_RES.X(2,:).^2);
    t = OPT_RES.t;
    Parameters = OPT_RES.param;
    hk = T / Parameters.N;
    cost_leg_e = 0;
    cost_leg_m = 0;
    cost_ankle_e = 0;
    cost_ankle_m = 0;
    R_leg = Parameters.R_leg;
    R_ankle = Parameters.R_ankle;
    mechA_leg = Parameters.mechA_leg;
    mechA_ankle = Parameters.mechA_ankle;
    for i = 1:Parameters.N - 1
        cost_leg_e = cost_leg_e + R_leg * u(1, i)^2 * hk;  %electrical
        cost_leg_m = cost_leg_m + max(0, u(1, i) * Parameters.transmission / mechA_leg * OPT_RES.X(6, i) * hk); %mechanical
        cost_ankle_e = cost_ankle_e + R_ankle * u(2, i)^2 * hk; %electrical
        cost_ankle_m = cost_ankle_m + max(0, u(2, i) * Parameters.transmission_ankle / (mechA_ankle * r(i)) * v(i) * hk); %mechanical
    end
    cost.leg_e = cost_leg_e;
    cost.leg_m = cost_leg_m;
    cost.ankle_e = cost_ankle_e;
    cost.ankle_m = cost_ankle_m;
    
    TorqueToStand = OPT_RES.param.m * OPT_RES.param.g / (OPT_RES.param.transmission/mechA_leg);
    cost.LegEtoStand = R_leg * TorqueToStand^2 * T;
    
    
    if plotifone == 1
        leg_e = R_leg * u(1,:).^2 * hk;
        leg_m = max(0, u(1,:) .* Parameters.transmission ./ mechA_leg .* OPT_RES.X(6, :) * hk);
        ankle_e = R_ankle * u(2,:).^2 * hk;
        ankle_m = max(0, u(2, :) .* Parameters.transmission_ankle / (mechA_ankle * r) * v * hk);
        
        figure;
        plot(t, leg_e); hold on;
        plot(t, leg_m); plot(t, ankle_e); plot(t, ankle_m);
        legend('Leg electrical', 'Leg mechanical', 'Ankle electrical', 'Ankle mechanical')
        
    end
end

