function [ cost ] = objectiveFunctionSymbolic( dv, p, sp)
%objective_function This function evaluates the objective function
%   This function uses Decision variables and calculates the objective
%   function or cost
    ind = 0:p.N-2;
    x = dv(ind*8+1);
    x(43);
    y = dv(ind*8+2);
%     r0 = dv(ind*8+3);
    dx = dv(ind*8+4);
    dy = dv(ind*8+5);
    dr0 = dv(ind*8+6);
    Tleg = dv(ind*8+7);
    Tankle = dv(ind*8+8);
    r = sqrt(x.^2 + y.^2);
    Tstance = dv(end-1);
    Tflight = dv(end);
    hk = [Tstance/p.Nstance * ones([1,p.Nstance]),...
        Tflight/p.Nflight * ones([1,p.Nflight])];
    cost_leg = 0;
    cost_ankle = 0;
    R_leg = sp(sp == 'R_leg');
    R_ankle = sp(sp == 'R_ankle');
    transmission = sp(sp == 'transmission');
    transmission_ankle = sp(sp == 'transmission_ankle');
%     R_leg = sp.R_leg;
%     R_ankle = sp.R_ankle;
%     transmission = sp.transmission;
%     transmission_ankle = sp.transmission_ankle;
    maxXzero = MikeMax(1);
    for i = 1:p.N-1
        cost_leg = cost_leg + ...
            R_leg * Tleg(i)^2 * hk(i) +... %electrical
            maxXzero(Tleg(i) * transmission * dr0(i) * hk(i)); %mechanical
        cost_ankle = cost_ankle + ...
            R_ankle * Tankle(i)^2 * hk(i) +... %electrical
            maxXzero(Tankle(i) * transmission_ankle *...
            (x(i) * dy(i) - y(i) * dx(i)) / (r(i)^2)  * hk(i)); %mechanical
    end
    cost = cost_leg + cost_ankle;
end