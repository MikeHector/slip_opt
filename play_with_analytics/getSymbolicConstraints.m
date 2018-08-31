function [ c, ceq, c_linear, A, b] = getSymbolicConstraints( dv, p, v )
%nonlcon This function determines the violation of the constraints
%   Using trapezoidal collocation, the violation of the constraints is
%   determined by plugging in the decision variables into formulated
%   dynamics such that we expect 0. The difference from that and 0 is our
%   violation

    %Stance collocation constraints
    ceq_stance = stanceConstraintsNew(dv, p, v);
    
    %Flight collocation constraints
    ceq_flight = flightConstraintsNew(dv, p, v);
    
    %Gather up vars
    A = zeros(1,size(dv,1));
    ind = 0:p.N-2;
    x = dv(ind*8+1);
    y = dv(ind*8+2);
    r0 = dv(ind*8+3);
    dx = dv(ind*8+4);
    dy = dv(ind*8+5);
    dr0 = dv(ind*8+6);
    Tleg = dv(ind*8+7);
    Tankle = dv(ind*8+8);
    r = sqrt(x.^2 + y.^2);

%     eqCon = zeros(9,3);  
    stanceEnd = p.Nstance;
    %Starting stance constraints
    %Fs starts at 0
    eqCon = r(1) - r0(1);
    %y - velocity + position. energy approach
    eqCon(end+1) = .5 * dy(1)^2 - v.g * (v.apex_height - y(1));
    %x - velocity
    c_linear = dx(1) - v.apex_velocity;
    A(1,getSymNum(dx(1))) = 1; b(1) = v.apex_velocity;
    
    
    %Ending stance constraints
    %Force in spring is 0
    eqCon(end+1) = r(stanceEnd) - r0(stanceEnd);
    %Liftoff energy
%     eqCon(5,1) = Parameters.g * (Parameters.apex_height + Parameters.deltah - y(stanceEnd)) - .5 * dy(stanceEnd)^2;
    %Apex velocity x
    c_linear(end+1) = dx(stanceEnd) - (v.apex_velocity);
    A(2,getSymNum(dx(stanceEnd))) = 1; b(2) = v.apex_velocity;

    %Transition Constraints
    %State at end of stance matches state at start of flight
%     eqCon(4:6,2) = dv(4:6,stanceEnd) - dv(4:6, stanceEnd+1);
    
    %Ending flight constraints - match initial stance translated in x
    qE = [y(end) r0(end) r0(end) dx(end) dy(end) dr0(end) Tleg(end) Tankle(end)];
    q1 = [y(1) r0(1) r0(1) dx(1) dy(1) dr0(1) Tleg(1) Tankle(1)];
    concatplease = qE - q1;
    c_linear = [c_linear , concatplease];
    for i = 1:7
        A(i+2,getSymNum(qE(i))) = 1; 
        A(i+2,getSymNum(q1(i))) = -1;
        b(i+2) = 0;
    end
        
%     %Lock the TD angle ~~~~~~~~~~~~~~~~~~~FIX~~~~~~~~~~~~~~~~~~~~~~~~
%     if ~isnan(Parameters.TD_disturb)
%         eqCon(9,end) = atan2(y(1),x(1)) - (Parameters.baseline_TDA + Parameters.TD_disturb);
%     end
    
    %Average velocity
    eqCon(end+1) = dx(1) - (x(end) - x(1)) / (dv(end-1) + dv(end));
    
    %Concatenate
    ceq = [ceq_stance; ceq_flight; eqCon'];
    c_linear = c_linear';

    %Inequality constraints
    %Ankle torque bounds
    ankle_bound = (v.lf .* y(1:stanceEnd) .* ...
                   v.k .* (r0(1:stanceEnd) - r(1:stanceEnd)))...
                   ./ (2 .* r(1:stanceEnd));
               
    Acon1 = Tankle(1:stanceEnd) * v.transmission_ankle - ankle_bound;
    Acon2 = -(Tankle(1:stanceEnd) * v.transmission_ankle + ankle_bound);
    
    c = [Acon1; Acon2];
    
    assert(size(A,1) == length(c_linear) -1, 'Number of linear constraints may have changed')
    
end
