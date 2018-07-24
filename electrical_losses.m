%MHector
%7.23.18
%Find Electrical losses parameter for leg and ankle actuators

%% Leg
maxPower_leg = 735;
noLoadSpeed_leg = 1300; %rpm
cycloidGearDrive = 16; %16:1 transmission
mechanicalGainSpeed = .25; %Ratio of leg joint velocity to toe velocity
mechanicalGainForce = 2; %Ratio of leg joint torque to toe GRF
uDotMaxLeg = noLoadSpeed_leg * 2 * pi/ 60 / cycloidGearDrive * mechanicalGainSpeed;
FmotorLeg = cycloidGearDrive * mechanicalGainForce;

R_leg = uDotMaxLeg^2 / maxPower_leg * FmotorLeg^2; %Multiply this by leg motor torque^2 to get electrical losses

%% Ankle
maxPower_ankle = 430;
noLoadSpeed_ankle = 2900; %rpm
harmonicGearDrive = 50; %50:1 transmission
lever = .07; %Lever from harmonic drive
uDotMaxAnkle = noLoadSpeed_ankle * 2 * pi/ 60/ harmonicGearDrive;
FmotorAnkle = harmonicGearDrive * lever;

R_ankle = uDotMaxAnkle^2 / maxPower_ankle * FmotorAnkle^2; %Multiply this by ankle motor torque^2 to get electrical losses



