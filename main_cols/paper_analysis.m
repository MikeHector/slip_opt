%MHector
%8.20.18
%Ankle questions

%% Test 1
%What is the total energy cost of adding 1 J of mechanical ankle energy
%throughout stance vs. 1 J of mechanical leg energy
clc; clear;
load('baseline5'); R = opt_results; clear opt_results;
hk = T / Parameters.N;
W_tot_leg = 1 + R.param.R_leg ./ ((R.Tleg .* R.dr0).^2 * hk);
