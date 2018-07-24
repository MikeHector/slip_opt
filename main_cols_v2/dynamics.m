function [ f_eval ] = dynamics( DecisionVars, Parameters, time )
%dynamics This function handles the discrete dynamics of the slip stance
%dynamics
%   This function takes a decision variable vector and returns the
%evaluated dynamics
    x = DecisionVars(1); y = DecisionVars(2); r0 = DecisionVars(3);
    dx = DecisionVars(4); dy = DecisionVars(5); dr0 = DecisionVars(6);
    Tleg = DecisionVars(7); Tankle = DecisionVars(8);
    
    if Parameters.disturbance_on == 1 && time > Parameters.disturbance_t_start && time < Parameters.disturbance_t_end
        Fp = Parameters.disturbance_f; %Perturbation 
    else
        Fp = 0;
    end
    
    %Common factors
    r = sqrt(x.^2 + y.^2);
    Fs = Parameters.k * (r0 - r);
    Ft = Tankle * Parameters.transmission_ankle / (Parameters.mechA_ankle * r) + Fp; %Force from ankle torque and perturbation
    Fd = Parameters.c * dr0;
    Fg = Fs - Tleg * Parameters.transmission / Parameters.mechA_leg - Fd;
    
    xdd = (Fs * x / r - Fd * x / r  + Ft * y / r) / Parameters.m;
    
    ydd = (Fs * y / r - Fd * y / r - Ft * x / r  - Parameters.m * Parameters.g) / Parameters.m;
    
    r0dd = -Fg / (Parameters.transmission^2 * Parameters.i_motor);
    
    f_eval = [dx, dy, dr0, xdd, ydd r0dd]';
end

