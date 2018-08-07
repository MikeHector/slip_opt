function [ c, ViolationMatrix ] = nonlinear_constraint_func( DecisionVars, Parameters )
%nonlcon This function determines the violation of the constraints
%   Using trapezoidal collocation, the violation of the constraints is
%   determined by plugging in the decision variables into formulated
%   dynamics such that we expect 0. The difference from that and 0 is our
%   violation

    ViolationMatrix = zeros(size(DecisionVars));
    T = DecisionVars(9,1);
    %Collocation constraints (dynamics)
    for i = 1:Parameters.N - 1
        time_i = T * (i - 1) / (Parameters.N - 1);
        time_i_plus_1 = T * i / (Parameters.N - 1);
        Fk = dynamics(DecisionVars(:, i), Parameters, time_i);
        Fk_plus_1 = dynamics(DecisionVars(:, i + 1), Parameters, time_i_plus_1);
        %Trapezoidal Collocation Constraints
        ViolationMatrix(:, i + 1) = ...
            [.5 * (time_i_plus_1 - time_i) * (Fk_plus_1 + Fk)...
            - DecisionVars(1:2*Parameters.dof, i + 1)...
            + DecisionVars(1:2*Parameters.dof, i); 0; 0; 0];
    end
    x = DecisionVars(1,:);
    y = DecisionVars(2,:);
    r0 = DecisionVars(3,:);
    dx = DecisionVars(4,:);
    dy = DecisionVars(5,:);
    dr0 = DecisionVars(6,:);
    Tleg = DecisionVars(7,:);
    Tankle = DecisionVars(8,:);
    r = sqrt(x.^2 + y.^2);
    
    %Starting constraints
    %r0 spring starts at nominal length
    ViolationMatrix(1,end + 1) = Parameters.r0_start - r0(1);
    %Fs starts at 0
    ViolationMatrix(2,end +1) = r(1) - r0(1);
    %Velocity of spring starts at 0
    ViolationMatrix(3,end) = dr0(1);
    %y - velocity + position. energy approach
    ViolationMatrix(4,end) = .5 * dy(1)^2 - Parameters.g * (Parameters.apex_height - y(1));
    %x - velocity
    ViolationMatrix(5,end) = dx(1) - Parameters.apex_velocity;

    %Ending Constraints
    %x - velocity end
    ViolationMatrix(6,end) = dx(end) - (Parameters.apex_velocity + Parameters.deltav);
    %Final apex height constraint
    ViolationMatrix(7,end) = Parameters.g * (Parameters.apex_height + Parameters.deltah - y(end)) - .5 * dy(end)^2;
    %End condition Fs =0 ~> r0 = r @ End of stance 
    ViolationMatrix(8,end) = r0(end) - r(end); 
    
    %Lock the TD angle
    if ~isnan(Parameters.TD_disturb)
        ViolationMatrix(9,end) = atan2(y(1),x(1)) - (Parameters.baseline_TDA + Parameters.TD_disturb);
    end
    
%     ViolationMatrix(1,end + 1) = y(1) - (y(end) - 0); %Dont need with new apex height constraint
%     ViolationMatrix(1,end+1) = x(1) + x(end); %Artificial AF

    %Inequality constraints
    %Ankle torque bounds
    ankle_bound = (Parameters.lf .* y .* ...
                   Parameters.k .* (r0 - r)) ./ (2 .* r);
    Acon1 = Tankle - ankle_bound;
    Acon2 = -(Tankle + ankle_bound);
    
    %/JUST FOR SEED FINDING!!!
%     shit = -x .* dy;
    %/END SEED FINDING
    
    A_vio = [Acon1; Acon2];
    
    c = [A_vio];
    
end
