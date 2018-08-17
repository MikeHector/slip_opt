%MHector
%8.16.18
%Find Electrical losses parameter for leg and ankle actuators
% v1 was so wrong RIP
clear; clc;
%Resistive power losses = T^2/km^2 -> Ploss = T^2 * R -> R = 1/km^2
%% Leg 
%Cassie leg motor: ILM 115x25
km_leg = 0.88; % Nm/sqrt(W)
R_leg = 1/km_leg^2;

%% Ankle
%Cassie ankle motor: ILM 85x13
km_ankle = 0.328; % Nm/sqrt(W)
R_ankle = 1/km_ankle^2;
