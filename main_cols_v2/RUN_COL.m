% MHector
% 4.26.18
% Collocation of slip through stance; v2 includes stance duration as a
% decision variable
function [optimized, opt_res] = RUN_COL(seed, leg_damping, apex_vel, apex_height, ankles_on, end_vel, disturbance_f, TD_angle) %TD_angle should be nan unless it is dictated
    %Define parameters
    close all
    param.m = 32; param.k = 4000; param.r0_start = .9; param.g= 9.81; 
    param.anklemax = 4.5; param.legmax = 12.2; param.N = 70;
    param.dof = 3; param.cntrl_dof = 2; param.i_motor = .0934/16^2; 
    param.transmission = 16; param.lf = .23; param.r0_min = .7; 
    param.r0_max = 1.3; param.timemax = 1; param.timemin = .2;
    param.c = leg_damping; param.fmincon_stuff = []; param.end_vel = end_vel;
    param.disturbance_on = 1; param.disturbance_f = disturbance_f;
    param.disturbance_t_start = 0.05; param.disturbance_t_end = .15;
    param.TD_angle = TD_angle; param.obj_func = 'Better approximation of energy';
    param.R_leg = 25.2; param.R_ankle = 214.5; param.apex_height_final = apex_height;
    param.mechA_leg = .2; param.mechA_ankle = .7; param.transmission_ankle = 50;
    
    assert(length(seed) == param.N,'Seed was not the expected dimension')
    %Define initial conditions
    param.apex_height = apex_height;
    param.apex_velocity = apex_vel;

    x0 = seed;

    %%%%%%%

    param.y_start = x0(2,1); param.x_start = x0(1,1); param.x_dot_start = x0(4,1); 
    param.y_dot_start = x0(5,1); %param.T = stance_duration;

    % plot(x0(1,:),x0(2,:))

    %Cost function
    cost_func = @(x) objective_function_2(x, param);

    %Set up nonlinear constraint (which we'll be calling a lot)
    nonlcon = @(x) nonlinear_constraint_func(x, param);

    %Set other constraints
    A = []; B = []; Aeq = []; Beq = [];

    %Set bounds
    ub = zeros(9, param.N) + inf; ub(3, :) = param.r0_max;
    ub(7, :) = param.legmax; ub(8, :) = param.anklemax; ub(9,1) = param.timemax;
    lb = zeros(9, param.N) - inf; lb(3, :) = param.r0_min;
    lb(7, :) = -param.legmax; lb(8, :) = -param.anklemax;
    lb(2,:) = zeros(1, param.N); lb(9,1) = param.timemin;

    if ankles_on == 0 %Turn off ankle torques
        ub(8,:) = 0; lb(8,:) = 0;
    end

    %Turn these on when you want to plot
%     plot_names.traj_fig = figure; plot_names.traj_plot = plot(1,1);
%     plot_names.cost_fig = figure; plot_names.cost_plot = plot(1,1);
%     plot_names.cost_plot.Marker = 'o';
%     plot_names.torque_fig = figure; plot_names.torque_plot_leg = plot(1,1);
%     plot_names.torque_plot_leg.Color = 'b';
%     figure(plot_names.torque_fig); hold on; plot_names.torque_plot_ankle = plot(1,1); hold off;
%     plot_names.torque_plot_ankle.Color = 'r';
%     outputfunc = @(x,optimValues,state) outfun(x,optimValues,state, plot_names);


    %Optimizer parameters
    param.options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1000000,...
                            'MaxIterations', 150000, 'ConstraintTolerance', 1e-3,...
                            'UseParallel', true, ...
                            'OptimalityTolerance',1e-6); %,'OutputFcn', outputfunc)


    %Optimize!
    [optimized, ~, param.flag, param.fmincon_stuff] = fmincon(cost_func,x0,A,B,Aeq,Beq,lb,ub,nonlcon,param.options);


    x = optimized(1,:); y = optimized(2,:); r0 = optimized(3,:);
    dx = optimized(4,:); dy = optimized(5,:); dr0 = optimized(6,:);
    Tleg = optimized(7,:); Tankle = optimized(8,:); T = optimized(9,1);
    t = linspace(0, T, param.N); r = sqrt(x.^2 + y.^2);

    opt_res = param;
    opt_res.x = x; opt_res.y = y; opt_res.r0 = r0; opt_res.dx = dx; opt_res.dy = dy;
    opt_res.dr0 = dr0; opt_res.Tleg = Tleg; opt_res.Tankle = Tankle; opt_res.T = T;
    opt_res.t = t; opt_res.r = r; opt_res.cost = cost_func(optimized);
    opt_res.X = optimized;

end