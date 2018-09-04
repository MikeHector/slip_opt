function [Fsx, Fsy] = getGRF(optres, oneForPlot)
    Nstance = optres.param.Nstance;
    Fs = optres.param.k * (optres.r0(1:Nstance) - optres.r(1:Nstance));
    xcomp = optres.x(1:Nstance)./optres.r(1:Nstance);
    ycomp = optres.y(1:Nstance)./optres.r(1:Nstance);
    Fsx = Fs .* xcomp;
    Fsy = Fs .* ycomp;
    if oneForPlot == 1
        figure;
        plot(optres.t(1:Nstance), Fsx); hold on;
        plot(optres.t(1:Nstance), Fsy);
        legend('X comp of GRF', 'Y comp of GRF')
        title('GRFs over time')
        xlabel('Time'); ylabel('Force')
        
    end
end