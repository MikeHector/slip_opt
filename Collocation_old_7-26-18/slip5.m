classdef slip5
    %SLIP class for traj opt
    %SLIP has leg motor
    %Includes ground damping
    %Has ankle
    
    properties
        m = 32; %mass
        k = 3000; %stiffness
        c = 0; %damping
        b = 100000000; %Ground damping
        lf = .2;
        thetac = 0 %Theta touchdown; \ <-- small pos angle
        r0_start = .9; %resting length of spring at start
        inertia_motor = .003;
        transmission = 50;
        g = 9.81; %Gravity
        
        ankle_torque_dot_limit = 400; %Nm/s
        leg_length_dot_limit = 20 %m/s;
        timestep = .005; %Time step for ode45
        step = 0;
        num_flights = 0;
        num_stances = 0;
        
        %Optimization stuff
        state_length = [];
        time_tape = []; %stance time span for torque tape
        leg_torque_tape = [];
        ankle_torque_tape = [];
        stance_end = [];
      
        x_body = 0; %initial x position at apex of flight
        dx_body = 1; %intial x velocity at apex of flight
        y_body = 1.1; %initial y position at apex of flight
        dy_body = 0; %initial y velocity at apex of flight
        r0 = .9;
        dr0 = 0;
        r_spring = .9;
        t = 0;
        r2 = 0;
        dr2 = 0;
        leg_angle = 0;
        
        kp = 30000;
        kd = 100000;
        
        dynamic_state = 0; %state of slip; 0 = flight, 1 = stance
        crash = 0; %0 = fine; 1 = crashing
        
        xtoe = 0;
        xtoe_new = [];
        cost = [];
        
        going_for_apex = 0;
        
        torque_track = [0,0];
        commanded_leg_torque = [];
        t_stance_start
        t_stance_end
        
        energy_consumed
    end
    
    methods
        
        function [opt_tape, flag] = optimize_leg(s, just_constraints)
            %This runs an optimization to find the optimal torque tape for
            %the leg length motor which satifies equilibrium conditions as
            %constraints and minimizes torque
            
            if just_constraints == 1
                cost_func = @(T) 0;    
                options = optimoptions('fmincon','Display','iter',...
                    'MaxFunctionEvaluations', 10000,'MaxIterations',5000,...
                    'ConstraintTolerance',1e-4,'UseParallel',true,...
                    'Algorithm','SQP','FiniteDifferenceStepSize', 1e-4,...
                    'FiniteDifferenceType', 'central');
            elseif just_constraints == 0
                cost_func = @(T) sum(T(2:end-1).^2) * 1e-5;
                options = optimoptions('fmincon','Display','iter',...
                    'MaxFunctionEvaluations', 10000,'MaxIterations',5000,...
                    'ConstraintTolerance',1e-2,'UseParallel',true,...
                    'Algorithm','SQP', 'FiniteDifferenceStepSize', 1e-4,...
                    'FiniteDifferenceType', 'central');
            end
            
            x0 = [s.thetac, s.leg_torque_tape, s.stance_end];
            A = [];
            B = [];
            Aeq = [];
            Beq = [];
            lb_thetac = -0.1;
            ub_thetac = 0.5;
            lb_leg = zeros(1, s.state_length) - 25;
            ub_leg = zeros(1, s.state_length) + 25;
            lb_t_stance_end = .4;
            ub_t_stance_end = .75;
            
            lb = [lb_thetac, lb_leg, lb_t_stance_end];
            ub = [ub_thetac, ub_leg, ub_t_stance_end];
            
            nonlcon = @(tape) s.find_nonlcon_leg(tape);            
            
            %Set up plots
            
            
            [opt_tape, ~, flag] = fmincon(cost_func,x0,A,B,Aeq,Beq,lb,ub,nonlcon,options);
        end
        
        function [c, ceq] = find_nonlcon_leg(s, tape)
            s.thetac = tape(1);
            stance_start = s.simSLIP(10,0,2,2,1,2).t_stance_start;
            s.leg_torque_tape = tape(2:s.state_length + 1);
            s.time_tape = linspace(stance_start, tape(end), s.state_length);
            if length(s.time_tape) ~= length(tape) -2
                disp('wtf2')
            end
            s_con = s.simSLIP(10,0,1,2,1,0);
            
            %Can it speed up
            ceq = [s_con.y_body(end) - s_con.y_body(1) ;...
                   s_con.dx_body(end) - s_con.dx_body(1) - 1];
               if isnan(ceq)
                   disp('wtf3')
               end
               
               
               
               if any(s_con.x_body > 1)
                   disp('didnt leave the ground')
               end
            c = 0;
            
            
            figure(11)
            hold on
            plot(s.time_tape, s.leg_torque_tape, 'r')
            title('Torque through stance')
            xlabel('Time')
            ylabel('Torque')
            
            figure(12)
            hold on
            plot(s_con.x_body, s_con.y_body)
            title('X vs Y of SLIP response')
            xlabel('X displacement')
            ylabel('Y displacement')
    
            drawnow
        end
        
        function [opt_tape, flag] = optimize_both(s, just_constraints)
            %This runs an optimization to find the optimal torque tape for
            %the leg length motor which satifies equilibrium conditions as
            %constraints and minimizes torque
            
            if just_constraints == 1
                cost_func = @(T) 0;    
                options = optimoptions('fmincon','Display','iter',...
                    'MaxFunctionEvaluations', 10000,'MaxIterations',5000,...
                    'ConstraintTolerance',1e-4,'UseParallel',true,...
                    'Algorithm','SQP','FiniteDifferenceStepSize', 1e-4,...
                    'FiniteDifferenceType', 'central');
            elseif just_constraints == 0
                cost_func = @(T) sum(T(2:end-1).^2) * 5e-6;
                options = optimoptions('fmincon','Display','iter',...
                    'MaxFunctionEvaluations', 10000,'MaxIterations',5000,...
                    'ConstraintTolerance',1e-2,'UseParallel',true,...
                    'Algorithm','SQP', 'FiniteDifferenceStepSize', 1e-4,...
                    'FiniteDifferenceType', 'central');
            end
            
            x0 = [s.thetac, s.leg_torque_tape, s.ankle_torque_tape, s.stance_end];
            A = [];
            B = [];
            Aeq = [];
            Beq = [];
            lb_thetac = -0.1;
            ub_thetac = 0.5;
            lb_leg = zeros(1, s.state_length) - 25;
            ub_leg = zeros(1, s.state_length) + 25;
            lb_ankle = zeros(1, s.state_length) - 10;
            ub_ankle = zeros(1, s.state_length) + 10;
            lb_t_stance_end = .4;
            ub_t_stance_end = .75;
            
            lb = [lb_thetac, lb_leg, lb_ankle, lb_t_stance_end];
            ub = [ub_thetac, ub_leg, ub_ankle, ub_t_stance_end];
            
            nonlcon = @(tape) s.find_nonlcon_both(tape);            
            
            [opt_tape, ~, flag] = fmincon(cost_func,x0,A,B,Aeq,Beq,lb,ub,nonlcon,options);
        end
        
        function [c, ceq] = find_nonlcon_both(s, tape)
            s.thetac = tape(1);
            stance_start = s.simSLIP(10,0,2,2,1,2).t_stance_start;
            s.leg_torque_tape = tape(2:s.state_length + 1);
            s.ankle_torque_tape = tape(s.state_length + 2: 2 * s.state_length + 1);
            s.time_tape = linspace(stance_start, tape(end), s.state_length);
            if length(s.time_tape) ~= length(tape)/2 -1
                disp('wtf2')
            end
            s_con = s.simSLIP(10,0,1,2,1,2);
            ceq = [s_con.y_body(end) - s_con.y_body(1);...
            s_con.dx_body(end) - s_con.dx_body(1);];
               if isnan(ceq)
                   disp('wtf3')
               end
            
            %Torque bound
            [~, r0_u_index] = unique(s_con.r0);
            r0_bound = interp1(s_con.t(r0_u_index), s_con.r0(r0_u_index), s_con.time_tape);
            [~, r_u_index] = unique(s_con.r_spring);
            r_bound = interp1(s_con.t(r_u_index), s_con.r_spring(r_u_index), s_con.time_tape);
            B = .47 * s_con.k * (r0_bound - r_bound) * s_con.lf / r_bound;
            
            

               
           c = [s_con.ankle_torque_tape - B;...
               -(B + s_con.ankle_torque_tape)];
               
               
               if any(s_con.x_body > 1)
                   disp('didnt leave the ground')
               end
            
            
            figure(11)
            hold on
            plot(s.time_tape, s.leg_torque_tape, 'r', s.time_tape, s.ankle_torque_tape, 'b')
            title('Torque through stance, red - leg, blue - ankle')
            xlabel('Time')
            ylabel('Torque')
            
            figure(12)
            hold on
            plot(s_con.x_body, s_con.y_body)
            title('X vs Y of SLIP response')
            xlabel('X displacement')
            ylabel('Y displacement')
    
            drawnow
        end
        
        function touchdown_angle_opt = find_EGB(s)
            %This runs an optimization to find thetac for classes
            %properties
            
            cost_func = @(theta) 0;
            
            x0 = s.thetac;
            A = [];
            B = [];
            Aeq = [];
            Beq = [];
            lb = -.1;
            ub = .5;

            nonlcon = @(theta) s.angle_opt_nonlincon(theta);
            
            options = optimoptions('fmincon','Display','iter',...
                'MaxFunctionEvaluations', 10000,'MaxIterations',5000,...
                'ConstraintTolerance',1e-3,'UseParallel',false,...
                'Algorithm','interior-point', 'FiniteDifferenceStepSize', 1e-4,...
                'FiniteDifferenceType', 'central');
            
            touchdown_angle_opt = fmincon(cost_func,x0,A,B,Aeq,Beq,lb,ub,nonlcon,options);
            
        end
        
        function cost = angle_opt_cost_func(s, theta_crit)
            s.thetac = theta_crit;
            disp(s.thetac)
            s_cost = simSLIP(s, 10, 0, 1, 0, 0, 0);
            cost = (s_cost.dx_body(1) - s_cost.dx_body(end))^2 +...
                   (s_cost.y_body(1) - s_cost.y_body(end))^2;

        end
        
        function [c, ceq] = angle_opt_nonlincon(s, theta_crit)
            s.thetac = theta_crit;
            s_con = simSLIP(s, 10, 0, 1, 0, 0, 0);
            ceq = [s_con.y_body(end) - s_con.y_body(1);...
                   s_con.dx_body(end) - s_con.dx_body(1)];
            c = 0;
        end
        
        function s = simSLIP(s, tspan, plot_results, one_bounce, leg_control, leg_motor, ankle_control)           
            while s.t < tspan(end) %Determine condition for full step maybe flight has a maximum not at the ends
                switch s.dynamic_state
                    case 0 %Flight
                        dyn = @(t,D) flightde(t,D,s);
                        if one_bounce == 1 && s.num_stances == 1
                            options = odeset('Events',@(t,D) flight_apex_event(t,D,s));
                            s.going_for_apex = 1;
                        else
                            options = odeset('Events',@(t,D) flight_event(t,D,s));
                        end
                        
                        if one_bounce == 2
                            options = odeset('Events',@(t,D) flight_event(t,D,s));
                            s.going_for_apex = 2;
                        end
                        q0 = [s.x_body(end)...
                              s.dx_body(end)...
                              s.y_body(end)...
                              s.dy_body(end)...
                              s.r0(end)...
                              0];
                          
                    case 1 %Stance
                        dyn = @(t,D) stancede(t,D,s, leg_control, leg_motor, ankle_control);
                        options = odeset('Events',@(t,D) stance_event(t,D,s),'RelTol',1e-8);
                        q0 = [s.x_body(end)...
                              s.dx_body(end)...
                              s.y_body(end)...
                              s.dy_body(end)...
                              s.r0(end)...
                              0];
                
                end
                tspan = s.t(end):s.timestep:tspan(end); %tspan for ode45
                [t_out,q_out,~,~,~] = ode45(dyn,tspan,q0,options); %Solve dyns
                s = track_data(s,t_out,q_out);
                s = determine_new_step(s, q_out);
                
%                 if plot_results == 1
%                     s = slip_plot(s);
%                 end
                
                s = check_for_crash(s);
                if s.crash == 1 
%                     disp('robot fell')
                    break
                end
                
                s = switchMode(s, plot_results, ankle_control); %Switch mode
                
                if s.going_for_apex == 1 
                    if plot_results == 1
                        disp('Made it to apex!')
                    end
                    %STOP THE PROGRAM
                    break
                end
                
                if s.going_for_apex == 2
                    s.t_stance_start = s.t(end);
                    break
                end
                
            end
        end

        function s = determine_new_step(s, q)
            %Find if a new step has occured
            if s.dynamic_state == 0
                max_index = find(q(:,3)==max(q(:,3)));
                if max_index ~= 1 && max_index ~= length(q(:,3))
                    %If max y index is in middle, slip reached an apex
                    %This is robust to case in which the leg length motor
                    %retracts faster than the robot falls, so it enters
                    %flight phase without really taking a step.
                    s.step = s.step + 1;
                end
            end
        end
        
        function s = switchMode(s, plot_results, ankle_control)
            switch s.dynamic_state
                case 0 %Flight
                    s.dynamic_state = 1; %Stance
                    if s.num_flights == 0
                        s.t_stance_start = s.t(end);
                    end
                    if plot_results == 1
                        disp('Starting stance')
                        disp(s.t(end))
                        s = slip_plot(s, ankle_control);
                    end
                    s.num_flights = s.num_flights + 1;
                case 1 %Stance
                    s.dynamic_state = 0; %Flight
                    s.t_stance_end = s.t(end);
                    if plot_results == 1
                        disp('Starting flight')
                        disp(s.t(end))
                        s = slip_plot(s, ankle_control);
                    end
                    s.num_stances = s.num_stances +1;
            end
        end
        
        function s = track_data(s, t, q)
            s.t = [s.t;t];
            s.x_body = [s.x_body; q(:,1)];
            s.dx_body = [s.dx_body; q(:,2)];
            s.y_body = [s.y_body; q(:,3)];
            s.dy_body = [s.dy_body; q(:,4)];
            s.r0 = [s.r0; q(:,5)];
            s.dr0 = [s.dr0; q(:,6)];
            
            switch s.dynamic_state
                
                case 0 %Data during flight
                s.r_spring = [s.r_spring; (zeros(size(q(:,1))) + s.r_spring(end))];
                s.xtoe = [s.xtoe; zeros(size(q(:,1))) + s.xtoe(end)];
                s.xtoe(end) = s.x_body(end) + s.r0(end) * sin(s.thetac);
                s.leg_angle = [s.leg_angle; (zeros(size(q(:,1))) + s.thetac)];
                
                case 1 %Data during stance
                xtoe_stance = zeros(size(q(:,1))) + s.xtoe(end);
                s.xtoe = [s.xtoe; xtoe_stance];
                s.r_spring = [s.r_spring; (sqrt((q(:,1) - s.xtoe(end)).^2 + q(:,3).^2) + 0)];
                s.leg_angle = [s.leg_angle; atan2(xtoe_stance - q(:,1), q(:,3))];
                
            end %switch
        end %track_data
        
        function s = check_for_crash(s)
            if any(s.y_body < .01)
                s.crash = 1;
            end
        end

        function [leg_motor_torque_at_t ] = leg_torque_finder(s, t)
        %This function is going to take in the time-scheduled tape of ankle
        %torques and current time, and output the linearly interpolated
        %leg torque at that time
            if length(s.time_tape) ~= length(unique(s.time_tape))
                leg_motor_torque_at_t = s.leg_torque_tape(1);
                disp('wtf')
            end
                
            if s.time_tape(end) < t
                leg_motor_torque_at_t = s.leg_torque_tape(end);
            else
                %This is the one that normally gets run
                leg_motor_torque_at_t = interp1(s.time_tape, s.leg_torque_tape, t);
            end
            
            if isnan(leg_motor_torque_at_t)
                leg_motor_torque_at_t = s.leg_torque_tape(1);
            end
        end %leg_torque_finder
            
        function [ankle_torque ] = ankle_torque_finder(s, t)
        %This function is going to take in the time-scheduled tape of ankle
        %torques and current time, and output the linearly interpolated
        %ankle torque at that time
            if length(s.time_tape) ~= length(unique(s.time_tape))
                ankle_torque = s.ankle_torque_tape(1);
                disp('wtf11')
            end

            if s.time_tape(end) < t
                ankle_torque = s.ankle_torque_tape(end);
%                 disp('wtf22')
            else
                %This is the one that normally gets run
                ankle_torque = interp1(s.time_tape, s.ankle_torque_tape, t);
            end
            
            if isnan(ankle_torque)
                ankle_torque = s.ankle_torque_tape(1);
                disp('wtf33')
            end
        end %ankle_torque_finder
        
        function leg_torque = leg_motor_PD(s, r_de, dr_de, t)
            %Outputs leg torque as a PD around the unsprung spring length,
            %length of the spring and change in length of the spring from 
            %the differential equation.
            leg_torque = s.kp * (r_de - s.r0_start) - s.kd * (dr_de);
            %Torque limits
            if leg_torque > 25
                leg_torque = 25;
            elseif leg_torque < -25
                leg_torque = -25;
            end
        end %leg_motor_PD
        
        function DE = flightde(t,D,s)
%             [leg_motor_torque] = torque_finder(t, s);
            DE = zeros(4,1);
            DE(1) = D(2);
            DE(2) = 0;
            DE(3) = D(4);
            DE(4) = -s.g;
            DE(5) = D(6);
            DE(6) = 0;
            %D(1)=x, D(2)=x', D(3)=y, D(4)=y', D(5) = r0, D(6) = r0'
        end %flightde
        
        function DE = stancede(t,D,s, leg_control, leg_motor, ankle_control)
            %Stance DE as a function of ankle torque tape       
            r_to_0 = sqrt((D(1) - s.xtoe(end)).^2 + (D(3)).^2);
            rdot = ((D(1) - s.xtoe(end)) .* D(2) + (D(3)) .* D(4)) ./ r_to_0;
            y_comp = (D(3)) / r_to_0;
            x_comp = (D(1) - s.xtoe(end)) / r_to_0;
            Fs = s.k * (D(5) - r_to_0);
            Fd = s.c * rdot;
            
            if ankle_control == 2 %Turn on ankle torque
                ankle_torque = s.ankle_torque_finder(t);
            elseif ankle_control == 0 %Turn off ankle torque
                ankle_torque = 0;
            elseif ankle_control == 1 %Some kind of controller
                disp('Please write a controller here')
            end
            
            Ft = ankle_torque ./ r_to_0;
            
            if leg_control == 1 %PD on r0
                leg_torque = s.leg_motor_PD(D(5), D(6), t);
                Fg = Fs - leg_torque * s.transmission - Fd;
            elseif leg_control == 2 %Tape
                leg_torque = s.leg_torque_finder(t);
                Fg = Fs - leg_torque * s.transmission - Fd;
            elseif leg_control == 0 %No torque
                Fg = (Fs - 0 * s.transmission - Fd);
            end

            DE = zeros(6,1);
            DE(1)=D(2);
            DE(2)=(Fs * x_comp + Ft * y_comp) / (s.m);
            DE(3)=D(4);
            DE(4)=(Fs * y_comp - Ft * x_comp - s.m * s.g) / (s.m);
            DE(5)=D(6);
            if leg_motor == 1
                DE(6)= -Fg / (s.transmission ^ 2 * s.inertia_motor);
            else
                DE(6) = 0;
            end
            %D(1)=x, D(2)=x', D(3)=y, D(4)=y', D(5)=r0, D(6)=r0'
        end %stancede2
    
        function [position,isterminal,direction] = stance_event(t,D,s)
            %Define the stance end event function
            spring_deflection = D(5) - sqrt((D(1) - s.xtoe(end)).^2 + (D(3)).^2);
            position = [spring_deflection; D(3)]; % The value that we want to be zero
            isterminal = [1; 1];  % Halt integration 
            direction = [0; 0];   % The zero can be approached from either direction
        end %Stance_event
    
        function [position,isterminal,direction] = flight_event(t,D,s)
            %Define the first flight end event
            yc = D(5) * cos(s.thetac);
            position = [D(3) - yc; D(3)]; % The value that we want to be zero
            isterminal = [1; 1];  % Halt integration 
            direction = [-1; 0];   % The zero can be approached from either direction
        end %flight_event
        
        function [position,isterminal,direction] = flight_apex_event(t,D,s)
            %Define the first flight end event
            position = [D(4); D(3)]; % The value that we want to be zero
            isterminal = [1; 1];  % Halt integration 
            direction = [0; 0];   % The zero can be approached from either direction
        end %flight_apex_event
        
        function s = animate(s, record_video)
            
            if record_video==1
                v=VideoWriter('vel_4','MPEG-4');
                v.FrameRate=10;
                open(v);
            end
            
            %Calculate important things
            r_to_0 = sqrt((s.x_body - s.xtoe).^2 + s.y_body.^2);
            Fs = s.k * (s.r0 - s.r_spring);
            t_stance_start_index = min(find(s.t == s.t_stance_start)) + 1;
            t_stance_end_index = max(find(s.t == s.t_stance_end));
            stance_points = t_stance_end_index - t_stance_start_index + 1;
            stance_time = linspace(s.t(t_stance_start_index), s.t(t_stance_end_index), stance_points)';
            ankle_torque_stance = interp1(s.time_tape, s.ankle_torque_tape, linspace(s.time_tape(1), s.time_tape(end), stance_points))';  
            ankle_torque = [zeros(t_stance_start_index,1); ankle_torque_stance; zeros(length(s.t) - t_stance_end_index - 1,1)];
            Fs_stance = Fs(t_stance_start_index:t_stance_end_index);
            r_to_0_stance = r_to_0(t_stance_start_index:t_stance_end_index);
            y_body_stance = s.y_body(t_stance_start_index:t_stance_end_index);  
            xcop_stance = (-ankle_torque_stance .* r_to_0_stance) ./ (Fs_stance .* y_body_stance);
            xcop = [zeros(t_stance_start_index, 1); xcop_stance; zeros(length(s.t) - t_stance_end_index - 1, 1)];
            leg_torque_stance = interp1(s.time_tape, s.leg_torque_tape, linspace(s.time_tape(1), s.time_tape(end), stance_points))'; 
            leg_torque = [zeros(t_stance_start_index, 1); leg_torque_stance; zeros(length(s.t) - t_stance_end_index - 1, 1)];
            
            for i = 1:length(xcop)
                if xcop(i) > .5 || xcop(i) < -.5
                    xcop(i) = 0;
                end
            end
            
            %X Y Animation
            slipbodyanimate=figure('Name','Slip body animation','pos',[1921,0,1921,1000]);
            figure(slipbodyanimate);
            ax1 = subplot(4,4,[1,5,9]);
            ax1.XLim = [-.1 2.5];
            ax1.YLim = [-.1 1.2];
            axis equal
            %COM circle
            theta = linspace(0, 2 * pi, 100);
            x_c = 0 + .1 * cos(theta);
            y_c = 0 + .1 * sin(theta);
            h(1) = patch(x_c, y_c, 'red');
            %Spring Line
            x_l = [-.01, .01, .01, -.01];
            y_l = [0, 0, -s.r_spring(1), -s.r_spring(1)];
            h(2) = patch(x_l, y_l, 'blue');
            tf = hgtransform('Parent', ax1);
            set(h, 'Parent', tf)
            rline = refline([0 0]);
            rline.Color = 'r';
            rline.LineStyle = '--';
            title('X Y SLIP')
            
            %COP Animation
            figure(slipbodyanimate)
            ax2 = subplot(4,4,2);
            ax2.XLim = [-.2 .2];
            ax2.YLim = [-.1 .2];
            x_foot = [-s.lf/2, s.lf/2, s.lf/2, -s.lf/2];
            y_foot = [0, 0, .05, .05];
            patch(x_foot, y_foot, 'blue') %foot
            arrowx = .25*[0, .1, 0, -.1];
            arrowy = .25*[0, -.2, -.1, -.2];
            arrow = patch(arrowx, arrowy, 'red');
            arrow_move = @(x_pos) arrowx + x_pos;
            title('COP')
            
            %COP plot
            figure(slipbodyanimate)
            ax3 = subplot(4,4,6);
            cop_plot = plot(0,0);
            ax3.XLim = [min(s.t), max(s.t)];
            ax3.YLim = [-s.lf/2 - .01, s.lf/2 + .01];
            cop_ref_upper = refline(ax3, [0, s.lf/2]); cop_ref_upper.Color = 'r'; cop_ref_upper.LineStyle = '--';
            cop_ref_lower = refline(ax3, [0, -s.lf/2]); cop_ref_lower.Color = 'r'; cop_ref_lower.LineStyle = '--';
            title('COP')
            
            %Ankle torque plot
            figure(slipbodyanimate)
            ax4 = subplot(4,4,[3,4,7,8]);
            ankle_torque_plot = plot(0,0);
            ax4.XLim = [min(s.t), max(s.t)];
            ax4.YLim = [min(ankle_torque), max(ankle_torque)];
            title('Ankle torque trajectory')
            
            %Leg animation
            figure(slipbodyanimate)
            ax5 = subplot(4,4,[10,14]);
            leg_anim = plot(0,0);
            ax5.XLim = [-.2, .2];
            ax5.YLim = [0, max(s.r0) + .1];
            %COM circle
            theta = linspace(0, 2 * pi, 100);
            x_c = 0 + .1 * cos(theta);
            y_c = 0 + .1 * sin(theta);
            circle = patch(x_c, y_c, 'r');
            new_circle_y = @(y) .1 * sin(theta) + y;
            %Spring
            x_l = [-.01, .01, .01, -.01];
            y_l = [0, 0, 1.4, 1.4];
            patch(x_l, y_l,'b')
            title('Spring length animation')
            
            %Leg torque plot
            figure(slipbodyanimate)
            ax6 = subplot(4,4, [11,12,15,16]);
            leg_torque_plot = plot(0,0);
            ax6.XLim = [min(s.t), max(s.t)];
            ax6.YLim = [min(leg_torque), max(leg_torque)];
            title('Leg torque trajectory')
            

            for q = 1:2
                for i = 1:length(s.x_body)
                    %Body animation stuff
                    %Repatch the spring line
                    h(2).Vertices(3:4,2) = -(s.r_spring(i) + 0);
                    %HG transform stuff
                    Tx = makehgtform('translate',[s.x_body(i),s.y_body(i),0],'zrotate',s.leg_angle(i));
                    set(tf, 'Matrix', Tx)
                    
                    %COP animation
                    %Repatch the arrow                       
                    arrow.Vertices(:,1) = arrow_move(xcop(i));
                    
                    %COP Plot
                    cop_plot.XData = s.t(1:i);
                    cop_plot.YData = xcop(1:i);
                    
                    %Ankle torque plot
                    ankle_torque_plot.XData = s.t(1:i);
                    ankle_torque_plot.YData = ankle_torque(1:i);
                    
                    %Leg length animation
                    circle.Vertices(:,2) = new_circle_y(s.r0(i));
                    
                    %Leg torque plot
                    leg_torque_plot.XData = s.t(1:i);
                    leg_torque_plot.YData = leg_torque(1:i);

                    drawnow
                    
                    if record_video==1
                        F=getframe(gcf);
                        writeVideo(v,F);
                    end
%                     pause(.005)
                end
            end
            
            if record_video == 1
                close(v)
            end
        
        end %animate
    
        
        function s = slip_plot(s, ankle_control)
            figure(74)
            plot(s.t, s.r0, 'r', s.t, s.r_spring,'b')
            title('t vs r0 in red, t vs r in blue')
            
            figure(23)
            plot(s.x_body,s.y_body, s.xtoe, 0, 'bo')
            title('x vs y of COM')
            axis equal
            
            figure(11)
            torquePD = s.leg_motor_PD(s.r0, s.dr0);
            plot(s.t, torquePD)
            title('Torque from PD')
            
%             figure(12)
%             plot(s.t, s.r2)
%             title('ground squish')

            figure(32)
            r_to_0 = sqrt((s.x_body - s.xtoe).^2 + s.y_body.^2);
            Fs = s.k * (s.r0 - s.r_spring);
            GRFx = Fs .* (s.x_body - s.xtoe) ./ r_to_0;
            GRFy = Fs .* (s.y_body) ./ r_to_0;
            plot(s.t, GRFx, 'g', s.t, GRFy, 'r')
            title('GRFs, green is GRFx; red is GRFy')
            
            if s.num_flights > 0 && ankle_control == 2
                t_stance_start_index = min(find(s.t == s.t_stance_start));
                t_stance_end_index = max(find(s.t == s.t_stance_end));
                stance_points = t_stance_end_index - t_stance_start_index + 1;
                stance_time = linspace(s.t(t_stance_start_index), s.t(t_stance_end_index), stance_points);
                ankle_torque = interp1(s.time_tape, s.ankle_torque_tape, linspace(s.time_tape(1), s.time_tape(end), stance_points));  
                figure(15)
                plot(stance_time, ankle_torque)
                title('Ankle torque over time')
                
                Fs_stance = Fs(t_stance_start_index:t_stance_end_index);
                r_to_0_stance = r_to_0(t_stance_start_index:t_stance_end_index);
                y_body_stance = s.y_body(t_stance_start_index:t_stance_end_index);
                xcop = (-ankle_torque' .* r_to_0_stance) ./ (Fs_stance .* y_body_stance);
                
                figure(17)
                plot(stance_time(3:end-3), xcop(3:end-3))
                axis([.2, .45, -.05, .05])
                title('Center of pressure through stance')
            end

            %Energy
            PEg = s.m * s.g* s.y_body;
            PEs = .5 * s.k * (s.r0 - s.r_spring).^2;
            KE = .5 * s.m * (s.dx_body.^2 + s.dy_body.^2);
%             NC = 0 * .5 * s.b * (diff(s.r2) ./ diff(s.t)).^2;
            NC = s.t(2:end) * 0;
            figure(86)
            plot(s.t, PEg, 'b', s.t, PEs, 'k', s.t, KE, 'g', s.t(2:end), PEg(2:end)+PEs(2:end)+KE(2:end)-NC, 'c',s.t(2:end), NC, 'm')
            title('Energy: Grav - blue, Spring - black, Kinetic - green, damper - magenta, Total - cyan')     
            
        end %animate

    end %methods
    
end

