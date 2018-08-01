function [out] = outfun(x,optimValues,state,plot_names)
%UNTITLED2 I want to plot stuff during fmincon in a not stupid way
%   Detailed explanation goes here
    out = false;

    if mod(optimValues.iteration,10) == 1 || optimValues.iteration < 30 || (optimValues.iteration > 500 && mod(optimValues.iteration, 100) == 1)
        %x is the current point at the iteration
        x_traj = x(1,:);
        y_traj = x(2,:);
        leg_torque = x(7,:);
        ankle_torque = x(8,:);
        
        %Trajectory plot
        
        plot_names.traj_plot.XData = x_traj;
        plot_names.traj_plot.YData = y_traj;
        
        %Cost plot
        if plot_names.cost_plot.YData == 1
            plot_names.cost_plot.XData = [];
            plot_names.cost_plot.YData = [];
        end
        plot_names.cost_plot.XData = [plot_names.cost_plot.XData, optimValues.iteration];
        plot_names.cost_plot.YData = [plot_names.cost_plot.YData, optimValues.fval];
        
        
        %Torque plot
        plot_names.torque_plot_leg.XData = x_traj;
        plot_names.torque_plot_leg.YData = leg_torque;
        
        plot_names.torque_plot_ankle.XData = x_traj;
        plot_names.torque_plot_ankle.YData = ankle_torque;
        
        drawnow
    end
end