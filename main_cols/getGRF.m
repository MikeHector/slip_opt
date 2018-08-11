function [Fsx, Fsy] = getGRF(optres, oneForPlot)
    Fs = optres.param.k * (optres.r0 - optres.r);
    xcomp = optres.x/optres.r;
    ycomp = optres.y/optres.r;
    Fsx = Fs * xcomp;
    Fsy = Fs * ycomp;
    if oneForPlot == 1
        figure;
        plot(optres.t, Fsx); hold on;
        plot(optres.t, Fsy);
        legend('X comp of GRF', 'Y comp of GRF')
    end
end