function [eqStance] = stanceConstraints(dv, p)
%stanceConstraints runs collocation constraints and others
%   Takes in full decision variables, parameters and outputs the inequality
%   constraints and equality constraints.

    eqStance = zeros(9,p.Nstance-1);
    Tstance = dv(9,1);
    %Collocation constraints (dynamics)
    for i = 1:p.Nstance-1
        time_i = Tstance * (i - 1) / (p.Nstance - 1);
        time_i_plus_1 = Tstance * i / (p.Nstance - 1);
        Fk = stanceDyn(dv(:, i), p);
        Fk_plus_1 = stanceDyn(dv(:, i + 1), p);
        %Trapezoidal Collocation Constraints
        eqStance(:, i) = ...
            [.5 * (time_i_plus_1 - time_i) * (Fk_plus_1 + Fk)...
            - dv(1:2*p.dof, i + 1)...
            + dv(1:2*p.dof, i); 0; 0; 0];
    end

end