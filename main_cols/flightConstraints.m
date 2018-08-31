function [eqFlight] = flightConstraints(dv, p)
%flightConstraints runs collocation constraints and others
%   Takes in full decision variables, parameters and outputs the inequality
%   constraints and equality constraints.

    eqFlight = zeros(9, p.Nflight-1);
    Tflight = dv(9,2);
    %Collocation constraints (dynamics)
    for i = p.Nstance:p.Nstance + p.Nflight - 2
        time_i = Tflight * (i - 1) / (p.Nflight - 1);
        time_i_plus_1 = Tflight * i / (p.Nflight - 1); %This isn't the exact time, but delta time is fine
        Fk = flightDyn(dv(:, i), p);
        Fk_plus_1 = flightDyn(dv(:, i + 1), p);
        %Trapezoidal Collocation Constraints
        eqFlight(:, i-p.Nstance+1) = ...%i+1?
            [.5 * (time_i_plus_1 - time_i) * (Fk_plus_1 + Fk)...
            - dv(1:2*p.dof, i + 1)...
            + dv(1:2*p.dof, i); 0; 0; 0];
    end
end

