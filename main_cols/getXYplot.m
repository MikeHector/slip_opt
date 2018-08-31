function [FullTraj] = getXYplot(OPTRES, plotifone)
    xStance = OPTRES.x;
    yStance = OPTRES.y;
    dxStance = OPTRES.dx;
    dyStance = OPTRES.dy;
    
    tTD = sqrt(2 * (yStance(1) - OPTRES.param.apex_height)/-OPTRES.param.g);
    xApex = xStance(1) - OPTRES.param.apex_velocity * tTD;
    
    tI = linspace(0,tTD, 40);
    xI = xApex + dxStance(1) * tI;
    yI = OPTRES.param.apex_height - .5 * OPTRES.param.g * tI.^2;
    dyI = -OPTRES.param.g * tI;
    
    tAF = max(roots([-.5 * OPTRES.param.g, dyStance(end), yStance(end) - (OPTRES.param.apex_height + OPTRES.param.deltah)]));
    xAF = xStance(end) + dxStance(end) * tAF;
    yAF = yStance(end) + dyStance(end) * tAF -.5 * OPTRES.param.g * tAF^2;
    
    tO = linspace(0, tAF, 40);
    xO = xStance(end) + dxStance(end) * tO;
    yO = yStance(end) + dyStance(end) * tO - .5 * OPTRES.param.g * tO.^2;
    
    if plotifone == 1
        figure; plot(xI,yI); hold on; plot(xStance, yStance); hold on; plot(xO, yO)
    end
    FullTraj.x = [xI, xStance, xO];
    FullTraj.y = [yI, yStance, yO];
    FullTraj.t = [tI, tI(end) + OPTRES.t, tI(end) + OPTRES.t(end) + tO];
    if plotifone == 1
        plot([xStance(1), xStance(end)],[yStance(1), yStance(end)], 'bo');
        title('XY trajectory including flight phase')
        xlabel('X Displacement'); ylabel('Y Displacement')
    end
    FullTraj.xStart = xStance(1);
    FullTraj.yStart = yStance(1);
    FullTraj.xEnd = xStance(end);
    FullTraj.yEnd = yStance(end);
    FullTraj.tTD = tTD;

end
    
    