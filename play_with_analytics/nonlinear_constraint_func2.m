function [ c, ceq, eqFlight] = nonlinear_constraint_func2( dv, Parameters )
%nonlcon This function determines the violation of the constraints
%   Using trapezoidal collocation, the violation of the constraints is
%   determined by plugging in the decision variables into formulated
%   dynamics such that we expect 0. The difference from that and 0 is our
%   violation

    %Stance collocation constraints
    eqStance = stanceConstraints(dv, Parameters);
    
    %Flight collocation constraints
    eqFlight = flightConstraints(dv, Parameters);
    
    %Names
    x = dv(1,:);
    y = dv(2,:);
    r0 = dv(3,:);
    dx = dv(4,:);
    dy = dv(5,:);
    dr0 = dv(6,:);
    Tleg = dv(7,:);
    Tankle = dv(8,:);
    r = sqrt(x.^2 + y.^2);
    
    eqCon = zeros(6,3);  
    stanceEnd = Parameters.Nstance;
    %Starting stance constraints
    %Fs starts at 0
    eqCon(1,1) = r(1) - r0(1);
    %y - velocity + position. energy approach
    eqCon(2,1) = .5 * dy(1)^2 - Parameters.g * (Parameters.apex_height - y(1));
    %x - velocity
    eqCon(3,1) = dx(1) - Parameters.apex_velocity;
    
    %Ending stance constraints
    %Force in spring is 0
    eqCon(4,1) = r(stanceEnd) - r0(stanceEnd);
    %Liftoff energy
%     eqCon(5,1) = Parameters.g * (Parameters.apex_height + Parameters.deltah - y(stanceEnd)) - .5 * dy(stanceEnd)^2;
    %Apex velocity x
    eqCon(5,1) = dx(stanceEnd) - (Parameters.apex_velocity);

    %Transition Constraints
    %State at end of stance matches state at start of flight
%     eqCon(4:6,2) = dv(4:6,stanceEnd) - dv(4:6, stanceEnd+1);
    
    %Ending flight constraints - match initial stance translated in x
    eqCon(6,1) = dv(2,end) - dv(2,1);
    eqCon(1:6,2) = dv(3:8,end) - (dv(3:8,1) +...
        [0 0 0 0 0 0]');
%     eqCon(1,3) = dv(8,end) - dv(8,1) + 0;
    
        
%     %Lock the TD angle ~~~~~~~~~~~~~~~~~~~FIX~~~~~~~~~~~~~~~~~~~~~~~~
%     if ~isnan(Parameters.TD_disturb)
%         eqCon(9,end) = atan2(y(1),x(1)) - (Parameters.baseline_TDA + Parameters.TD_disturb);
%     end
    
    %Average velocity
    eqCon(1,3) = dx(1) - (x(end) - x(1)) / (dv(9,1) + dv(9,2));
    
    %Concatenate
    ceq = [eqStance, eqFlight, eqCon];
    
    %Inequality constraints
    %Ankle torque bounds
    ankle_bound = (Parameters.lf .* y(1:stanceEnd) .* ...
                   Parameters.k .* (r0(1:stanceEnd) - r(1:stanceEnd)))...
                   ./ (2 .* r(1:stanceEnd));
               
    Acon1 = Tankle(1:stanceEnd) * Parameters.transmission_ankle - ankle_bound;
    Acon2 = -(Tankle(1:stanceEnd) * Parameters.transmission_ankle + ankle_bound);
    
    c = [Acon1; Acon2];
    
end
