%MHector
%7.25.18
%Coll sturcture builder

function [ColStrucArray] = ColStrucBuilderTest()
%Damping
colStruc.direction = {'up'};
colStruc.varName = 'c';
colStruc.deltaVar = 5;
colStruc.varMax = 2500;
colStruc.varMin = -1;
colStruc.var = 150; %Initial variable value
ColStrucArray.Damping = colStruc;
clear colStruc

%Changing Apex Velocity
colStruc.direction = {'up','down'};
colStruc.varName = 'apex_velocity';
colStruc.deltaVar = .01;
colStruc.varMax = .01;
colStruc.varMin = .01;
colStruc.var = 0;
ColStrucArray.ApexVelocity = colStruc;
clear colStruc

%Changing Force Disturbance
colStruc.direction = {'up', 'down'};
colStruc.varName = 'disturbance_f';
colStruc.deltaVar = .5;
colStruc.varMax = .5;
colStruc.varMin = -.5;
colStruc.var = 0;
ColStrucArray.ForceDisturbance = colStruc;
clear colStruc

%Changing Touch Down Angle Error
colStruc.direction = {'up','down'};
colStruc.varName = 'TD_disturb';
colStruc.deltaVar = .005;
colStruc.varMax = .005;
colStruc.varMin = -.005;
colStruc.var = 0;
ColStrucArray.TouchDownAngleError = colStruc;
clear colStruc

%Changing Velocity Between Apexes
colStruc.direction = {'up','down'};
colStruc.varName = 'deltav';
colStruc.deltaVar = .01;
colStruc.varMax = .01;
colStruc.varMin = -.01;
colStruc.var = 0;
ColStrucArray.VelocityBetweenApexes = colStruc;
clear colStruc

%Changing Heights Between Apexes
colStruc.direction = {'up','down'};
colStruc.varName = 'deltah';
colStruc.deltaVar = .005;
colStruc.varMax = .005;
colStruc.varMin = .005;
colStruc.var = 0;
ColStrucArray.VelocityBetweenApexes = colStruc;
clear colStruc

%Changing motor inertias
colStruc.direction = {'up'};
colStruc.varName = 'i_motor';
colStruc.deltaVar = (.095-.003)/12;
colStruc.varMax = .01;
colStruc.varMin = .002;
colStruc.var = .003;
ColStrucArray.DifferentMotorInertias = colStruc;
clear colStruc

end
