function [c] = stanceConstraintsNew(dv, p, sp)
%stanceConstraints runs collocation constraints and others
%   Takes in full decision variables, parameters and outputs the inequality
%   constraints and equality constraints.

%     eqStance = zeros(9,p.Nstance-1);
    c = [];
    %Collocation constraints (dynamics)
    for i = 1:p.Nstance-1
        time_i = dv(end-1) * (i - 1) / (p.Nstance - 1);
        time_i_plus_1 = dv(end-1) * i / (p.Nstance - 1);
        Fk = stanceDynNew(dv(8*(i-1)+1:8*(i)), sp);
        Fk_plus_1 = stanceDynNew(dv(8*(i)+1:8*(i+1)), sp);
        %Trapezoidal Collocation Constraints
        c = [c;...
            .5 * (time_i_plus_1 - time_i) * (Fk_plus_1 + Fk)...
            - dv(8*(i)+1:8*(i+1)-2)...
            + dv(8*(i-1)+1:8*i-2)];
    end

end