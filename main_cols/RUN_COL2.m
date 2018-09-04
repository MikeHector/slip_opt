% MHector
% 4.26.18
% Collocation of slip through stance; v2 includes stance duration as a
% decision variable
function [optimized, opt_res] = RUN_COL2(seed, param, smooth) %TD_angle should be nan unless it is dictated
    assert(length(seed) == param.N-1,'Seed was not the expected dimension')
    assert(param.Nstance+param.Nflight == param.N, 'Collocation points of stance and flight do not add up to param.N!')
    
    %Cost function
    cost_func = @(x) objective_function_3(x, param, smooth);

    %Set up nonlinear constraint (which we'll be calling a lot)
    nonlcon = @(x) nonlinear_constraint_func2(x, param);

    %Set other constraints
    A = []; B = []; Aeq = []; Beq = [];

    %Set bounds
    ub = zeros(9, param.N-1) + inf; %Fill out ub matrix
    ub(3, :) = param.r0_max; %Leg length max
    ub(7, :) = param.legmax; %Leg torque max
    ub(8, :) = param.anklemax; %Ankle torque max
    ub(9,1) = param.timemax; %Stance time max
    ub(9,2) = param.timemax; %Flight time max
    ub(9,3) = 4;
    lb = zeros(9, param.N-1) - inf; %Fill out lb matrix
    lb(3, :) = param.r0_min; %Leg length min
    lb(7, :) = -param.legmax; %Leg torque min
    lb(8, :) = -param.anklemax; %Ankle torque min
    lb(2,:) = zeros(1, param.N-1); %Y minimum of 0
    lb(9,1) = param.timemin; %Stance time min
    lb(9,2) = param.timemin; %Flight time min
    lb(9,3) = .5;

    if param.ankles_on == 0 %Turn off ankle torques
        ub(8,:) = 0; lb(8,:) = 0;
    end

    %Optimizer parameters
    param.options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1000000,...
                            'MaxIterations',15000000, 'ConstraintTolerance', 1e-3,...
                            'UseParallel', true, ...
                            'OptimalityTolerance',1e-3); %,'OutputFcn', outputfunc)

    %Optimize!
    [optimized, ~, param.flag, param.fmincon_stuff] = fmincon(cost_func,seed,A,B,Aeq,Beq,lb,ub,nonlcon,param.options);

    x = optimized(1,:); y = optimized(2,:); r0 = optimized(3,:);
    dx = optimized(4,:); dy = optimized(5,:); dr0 = optimized(6,:);
    Tleg = optimized(7,:); Tankle = optimized(8,:); Tstance = optimized(9,1);
    Tflight = optimized(9,2); apex_height = optimized(9,3); apex_velocity = optimized(4,1);
    tstance = linspace(0, Tstance, param.Nstance); average_velocity = param.average_velocity;
    tflight =  linspace(Tstance, Tstance+Tflight, param.Nflight);
    t = [tstance, tflight(2:end)];
    r = sqrt(x.^2 + y.^2); 

    opt_res.param = param;
    opt_res.x = x; opt_res.y = y; opt_res.r0 = r0; opt_res.dx = dx; opt_res.dy = dy;
    opt_res.dr0 = dr0; opt_res.Tleg = Tleg; opt_res.Tankle = Tankle; opt_res.Tstance = Tstance;
    opt_res.t = t; opt_res.r = r; opt_res.cost = cost_func(optimized);
    opt_res.X = optimized; opt_res.Tflight = Tflight; opt_res.param.apex_height = apex_height;
    opt_res.param.apex_velocity = apex_velocity; opt_res.param.average_velocity = average_velocity;

end