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
    %r0 spring starts undeflected; change in r0 is 0
    %Column 21
    ViolationMatrix(3,end + 1) = Parameters.r0_start - r0(1);
    ViolationMatrix(1,end) = Parameters.r0_start - r(1);
    ViolationMatrix(6,end) = dr0(1);

    %y - velocity + position. energy approach
    ViolationMatrix(2,end +1) = .5 * dy(1)^2 - Parameters.g * (Parameters.apex_height - y(1));
    
    %Lock the TD angle
    if ~isnan(Parameters.TD_angle)
        ViolationMatrix(2,end +1) = atan2(y(1),x(1)) - Parameters.TD_angle;
    end
    
    %x - velocity
    %Column 24
    ViolationMatrix(4,end + 1) = dx(1) - Parameters.apex_velocity;
    
    %x - position (with y-position and r0, would this overdefine the sys?)
    %This would give a more symmetric slip resposne
%     ViolationMatrix(1,end+1) = x(1) + x(end);
    
    %Symmetry constraint
%     if mod(Parameters.N,2) ~= 1
%         disp('NEED AN ODD NUMBER OF COLLOCATION POINTS')
%     end
%     for i = 1:((Parameters.N-1)/2)
%         symm_vio_x(i) = x(i) + x(Parameters.N - (i - 1));
%         symm_vio_y(i) = y(i) - y(Parameters.N - (i - 1));
%     end
%     ViolationMatrix(1,end:end + length(symm_vio_x) -1) = symm_vio_x;
%     ViolationMatrix(2,end - length(symm_vio_x) + 1:end) = symm_vio_y;
    
% %     If we want the TD angle to be EGB angle
%     ViolationMatrix(1,end+1) = Parameters.EGB - (atan2(y(1), x(1)) - pi/2);
    
%     %Artificial? dy at center of stance is 0
%     ViolationMatrix(5, end+1) = Parameters
    
    % Equilibrium gait constraints
    %Velocity EQ gait constraint
    %Column 25
    ViolationMatrix(4,end + 1) = dx(end) - Parameters.end_vel; %xvel - same
    ViolationMatrix(5,end + 1) = dy(1) + dy(end); %yvel - change sign
    
    %Height Eq gait constraint
    %Column 26
    ViolationMatrix(2,end + 1) = y(1) - (y(end) - 0);
    
    %End condition Fs =0 ~> r0 = r @ End of stance 
    %Column 27
    ViolationMatrix(3,end + 1) = r0(end) - r(end); 
    
    %Inequality constraints
    %Ankle torque bounds
    ankle_bound = (Parameters.lf .* y .* ...
                   Parameters.k .* (r0 - r)) ./ (2 .* r);
    Acon1 = Tankle - ankle_bound;
    Acon2 = -(Tankle + ankle_bound);
    A_vio = [Acon1; Acon2];
    
    c = [A_vio];

end
