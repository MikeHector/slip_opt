%MHector
%7.17.18
clc; clear;

n = 50;
m = 1; %percent

%Damping
directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\damping_results\opt*';
variable_name = 'c';

%Velocity
directory = 'C:\\Users\mike-\Documents\DRL\collocation\opt_results\velocity_results\vel*';
variable_name = 'apex_velocity';



[good, needs_improvement] = space_func(directory, variable_name, n, m);

needs_improvement


