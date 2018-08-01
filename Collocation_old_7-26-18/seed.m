%MHector 5.16.18
%Make a slip seed for collocation
%Take in stance time, (# of collocation points), and parameters, then
%spits out stance trajectories - x,y,dx,dy,r0,dr0,Tl,Ta interpolated by
%time at collocation points
clear; clc

% inputs
N = 20;
param.m = 32; param.k = 3000; param.r0_start = .9; param.g= 9.81; 
param.anklemax = 20; param.legmax = 25; param.T =2; param.N = 20;
param.dof = 3;param.i_motor = .003; param.transmission = 50; 
param.y_start = 1.1; param.x_start = -.2; param.x_dot_start = 1; 
param.y_dot_start = -.2; param.r0_max = 1.3; param.r0_min = .7; 
param.lf = .23;

% func
s = slip5;

s.m = param.m; s.k = param.k; s.r0_start = param.r0_start;
s.inertia_motor = param.i_motor; s.transmission = param.transmission;

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
stance_Ta = zeros(size(q.t));

start.t = linspace(stance_t(1),stance_t(end),N);
start.T1 = interp1(stance_t, stance_Tl, start.t);
start.x = interp1(stance_t, stance_x, start.t);
start.y = interp1(stance_t, stance_y, start.t);
start.dx = interp1(stance_t, stance_dx, start.t);
start.dy = interp1(stance_t, stance_dy, start.t);
start.r0 = interp1(stance_t, stance_r0, start.t);
start.dr0 = interp1(stance_t, stance_dr0, start.t);
start.Ta = interp1(stance_t, stance_Ta, start.t);









