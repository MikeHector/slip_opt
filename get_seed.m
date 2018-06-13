function [seed] = get_seed(param)
%get_seed Takes number of collocation parameters and generates a seed
%   Takes number of collocation points and parameters and runs a slip
%with EGB generated thetac and stiff PD on leg length to get a seed of
%decision variables 

s = slip5;

s.m = param.m; s.k = param.k; s.r0_start = param.r0_start;
s.r0 = param.r0_start; s.lf = param.lf;
s.inertia_motor = param.i_motor; s.transmission = param.transmission;
s.dx_body = param.apex_velocity; s.y_body = param.apex_height;

s.thetac = s.find_EGB;
q = s.simSLIP(10,0,1,1,1,0);

%Calculate leg torque
leg_torque = q.leg_motor_PD(q.r0, q.dr0, q.t);

%Get stance stuff
stance_start_index = max((find(q.t == q.t_stance_start)));
stance_end_index = min((find(q.t == q.t_stance_end)));
stance_Tl = leg_torque(stance_start_index:stance_end_index);
stance_t = q.t(stance_start_index:stance_end_index)- q.t_stance_start;
stance_x = q.x_body(stance_start_index:stance_end_index);
stance_y = q.y_body(stance_start_index:stance_end_index);
stance_dx = q.dx_body(stance_start_index:stance_end_index);
stance_dy = q.dy_body(stance_start_index:stance_end_index);
stance_r0 = q.r0(stance_start_index:stance_end_index);
stance_dr0 = q.dr0(stance_start_index:stance_end_index);
stance_Ta = zeros(size(stance_t));

x_offset = stance_x(stance_y == min(stance_y));

seed.t = linspace(stance_t(1),stance_t(end),param.N);
seed.Tl = interp1(stance_t, stance_Tl, seed.t);
seed.x = interp1(stance_t, stance_x, seed.t) - x_offset;
seed.y = interp1(stance_t, stance_y, seed.t);
seed.dx = interp1(stance_t, stance_dx, seed.t);
seed.dy = interp1(stance_t, stance_dy, seed.t);
seed.r0 = interp1(stance_t, stance_r0, seed.t);
seed.dr0 = interp1(stance_t, stance_dr0, seed.t);
seed.Ta = interp1(stance_t, stance_Ta, seed.t);
seed.stance_dur = [q.t_stance_end - q.t_stance_start, zeros(1, param.N-1)];
seed.EGB = s.thetac;
end

