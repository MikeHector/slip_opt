function [c] = flightConstraintsNew(dv, p, sp)
%flightConstraints runs collocation constraints and others
%   Takes in full decision variables, parameters and outputs the inequality
%   constraints and equality constraints.

%     eqFlight = zeros(9, p.Nflight-1);
%     Tflight = dv(9,2);
    c = [];
    %Collocation constraints (dynamics)
    for i = p.Nstance:p.Nstance + p.Nflight - 2
        time_i = dv(end) * (i - 1) / (p.Nflight - 1);
        time_i_plus_1 = dv(end) * i / (p.Nflight - 1); %This isn't the exact time, but delta time is fine
        Fk = flightDynNew(dv(8*(i-1)+1:8*(i)), sp);
        Fk_plus_1 = stanceDynNew(dv(8*(i)+1:8*(i+1)), sp);
        %Trapezoidal Collocation Constraints
        c = [c;...
            .5 * (time_i_plus_1 - time_i) * (Fk_plus_1 + Fk)...
            - dv(8*(i)+1:8*(i+1)-2)...
            + dv(8*(i-1)+1:8*i-2)];
    end
end

