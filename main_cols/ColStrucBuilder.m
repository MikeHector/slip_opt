%MHector
%7.25.18
%Coll sturcture builder

function [ColStrucArray] = ColStrucBuilder()
%Damping
colStruc.direction = {'up'};
colStruc.varName = 'c';
colStruc.deltaVar = 1;
colStruc.varMax = 1000;
colStruc.varMin = -1;
colStruc.var = 0; %Initial variable value
ColStrucArray.Damping = colStruc;
clear colStruc

%Changing Apex Velocity
colStruc.direction = {'up','down'};
colStruc.varName = 'apex_velocity';
colStruc.deltaVar = .01;
colStruc.varMax = 5;
colStruc.varMin = .1;
colStruc.var = 0;
ColStrucArray.ApexVelocity = colStruc;
clear colStruc

%Changing Force Disturbance
colStruc.direction = {'up'};
colStruc.varName = 'disturbance_f';
colStruc.deltaVar = .5;
colStruc.varMax = 50;
colStruc.varMin = -1;
colStruc.var = 0;
ColStrucArray.ForceDisturbance = colStruc;
clear colStruc

%Changing Touch Down Angle Error
colStruc.direction = {'up','down'};
colStruc.varName = 'TD_disturb';
colStruc.deltaVar = .005;
colStruc.varMax = 1;
colStruc.varMin = -1;
colStruc.var = 0;
ColStrucArray.TouchDownAngleError = colStruc;
clear colStruc

%Changing Velocity Between Apexes
colStruc.direction = {'up','down'};
colStruc.varName = 'deltav';
colStruc.deltaVar = .01;
colStruc.varMax = 3;
colStruc.varMin = -.75;
colStruc.var = 0;
ColStrucArray.VelocityBetweenApexes = colStruc;
clear colStruc

%Changing Heights Between Apexes
colStruc.direction = {'up','down'};
colStruc.varName = 'deltah';
colStruc.deltaVar = .005;
colStruc.varMax = .5;
colStruc.varMin = -.1;
colStruc.var = 0;
ColStrucArray.VelocityBetweenApexes = colStruc;
clear colStruc

end
