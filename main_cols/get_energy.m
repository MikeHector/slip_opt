function [cost] = get_energy(OPT_RES, plotifone)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    u = OPT_RES.X(7:8, :);
    T = OPT_RES.X(9,1);
    v = sqrt(OPT_RES.X(4,:).^2 + OPT_RES.X(5,:).^2);
    r = sqrt(OPT_RES.X(1,:).^2 + OPT_RES.X(2,:).^2);
    Parameters = OPT_RES.param;
    hk = T / Parameters.N;
    cost_leg_e = 0;
    cost_leg_m = 0;
    cost_ankle_e = 0;
    cost_ankle_m = 0;
    R_leg = Parameters.R_leg;
    R_ankle = Parameters.R_ankle;
    mechA_leg = 1;
    mechA_ankle = 1;
    for i = 1:Parameters.N
        cost_leg_e = cost_leg_e + R_leg * u(1, i)^2 * hk;  %electrical
        cle(i) = R_leg * u(1, i)^2;  %electrical
        cost_leg_m = cost_leg_m + max(0, u(1, i) * Parameters.transmission / mechA_leg * OPT_RES.X(6, i) * hk); %mechanical
        clm(i) = max(0, u(1, i) * Parameters.transmission / mechA_leg * OPT_RES.X(6, i)); %mechanical
        cost_ankle_e = cost_ankle_e + R_ankle * u(2, i)^2 * hk; %electrical
        cae(i) = R_ankle * u(2, i)^2; %electrical
        cost_ankle_m = cost_ankle_m + max(0, u(2, i) * Parameters.transmission_ankle *...
        (OPT_RES.x(i) * OPT_RES.dy(i) -OPT_RES.y(i) * OPT_RES.dx(i))/r(i)^2 * hk); %mechanical
        cam(i) = max(0, u(2, i) * Parameters.transmission_ankle *...
        (OPT_RES.x(i) * OPT_RES.dy(i) -OPT_RES.y(i) * OPT_RES.dx(i))/r(i)^2); %mechanical
    end
    cost.leg_e = cost_leg_e;
    cost.leg_m = cost_leg_m;
    cost.ankle_e = cost_ankle_e;
    cost.ankle_m = cost_ankle_m;
    
    TorqueToStand = OPT_RES.param.m * OPT_RES.param.g / (OPT_RES.param.transmission/mechA_leg);
    cost.LegEtoStand = R_leg * TorqueToStand^2 * T;
    
    if plotifone == 1
        figure;
        plot(OPT_RES.t/OPT_RES.t(end), cle); hold on;
        plot(OPT_RES.t/OPT_RES.t(end), clm); plot(OPT_RES.t/OPT_RES.t(end), cae); plot(OPT_RES.t/OPT_RES.t(end), cam);
        stand = refline(0,cost.LegEtoStand); stand.LineStyle = '--';
        legend('Leg electrical', 'Leg mechanical', 'Ankle electrical', 'Ankle mechanical', 'Power to stand')
        title('Actuator power contributing to energy consumption')
        xlabel('Time'); ylabel('Actuator power')
    end
end

