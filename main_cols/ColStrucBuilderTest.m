%MHector
%7.25.18
%Coll sturcture builder

function [ColStrucArray] = ColStrucBuilderTest()
%Damping
colStruc.direction = {'up'};
colStruc.varName = 'c';
colStruc.deltaVar = 50;
colStruc.varMax = 51;
colStruc.varMin = -1;
colStruc.var = 0; %Initial variable value
ColStrucArray.Damping = colStruc;
clear colStruc

%Changing Apex Velocity
colStruc.direction = {'up','down'};
colStruc.varName = 'apex_velocity';
colStruc.deltaVar = .5;
colStruc.varMax = 1.6;
colStruc.varMin = .4;
colStruc.var = 1;
ColStrucArray.ApexVelocity = colStruc;
clear colStruc

%Changing Force Disturbance
colStruc.direction = {'up', 'down'};
colStruc.varName = 'disturbance_f';
colStruc.deltaVar = 20;
colStruc.varMax = 21;
colStruc.varMin = -21;
colStruc.var = 0;
ColStrucArray.ForceDisturbance = colStruc;
clear colStruc

%Changing Touch Down Angle Error
colStruc.direction = {'up','down'};
colStruc.varName = 'TD_disturb';
colStruc.deltaVar = .05;
colStruc.varMax = .06;
colStruc.varMin = -.06;
colStruc.var = 0;
ColStrucArray.TouchDownAngleError = colStruc;
clear colStruc

%Changing Velocity Between Apexes
colStruc.direction = {'up','down'};
colStruc.varName = 'deltav';
colStruc.deltaVar = .5;
colStruc.varMax = .6;
colStruc.varMin = -.6;
colStruc.var = 0;
ColStrucArray.VelocityBetweenApexes = colStruc;
clear colStruc

%Changing Heights Between Apexes
colStruc.direction = {'up','down'};
colStruc.varName = 'deltah';
colStruc.deltaVar = .1;
colStruc.varMax = .11;
colStruc.varMin = -.11;
colStruc.var = 0;
ColStrucArray.HeightBetweenApexes = colStruc;
clear colStruc

% %Changing motor inertias
% colStruc.direction = {'up'};
% colStruc.varName = 'i_motor';
% colStruc.deltaVar = (.095-.003)/12;
% colStruc.varMax = .01;
% colStruc.varMin = .002;
% colStruc.var = .003;
% ColStrucArray.DifferentMotorInertias = colStruc;
% clear colStruc

end
