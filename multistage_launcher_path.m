close all
clc


%%%%%
%%
%% This script aim is finding following itemes
%% for VA_Ms_practicas_ejecicio1_22-23_enunciado_V4
%% 2b) Show that we are achieving results under error tolerance 
%%     during two orbits around earth.
%%     (error tolerance: 5% for major axis and excentricity less than 0.03)
%% 2c) Excess of propellant (if any)
%% 2d) Plot gravity, guiding and drag speed lost against ideal speed. 
%%
%%%%%
%Parameters
params = paramSet('T', [2718 798 83 0]*1e3,...
				 'Mp', [167.5 34.6 10.8 0]*1e3,...
				 'M', [185.015 38 12 3.38]*1e3,...
				 'tb',[150 125 600 inf],...
				 'Sref',pi*(5)^2/4,...
				 'g0',9.81);
Rt = 6378e3 %earth radius m 
earth_gravitacional_constant = 9.81*Rt^2; %m³/s²
startTimeGT = 45; %s time when gravity turn starts
endTimeGT = 700; %s time when gravity turn 'ends' 
startTimeGuiding = endTimeGT; % time when guiding by nozzle law enters 
target_height = 650e3; %m

%speed needed for a successful circle orbit around earth at 1230 km.
orbital_speed = sqrt(earth_gravitacional_constant/(Rt + target_height)); %m/s


%Gamma shooting, solution 3.08 degrees. See
%multistage_launcher_gravity_turn_angle.m
%% Apartado 2             
vertical_path = intMultiStage2DWithDrag([0 startTimeGT],[1e-7; 0; 1e-7; 0], params);
    params = paramSet('T', [2718 798 83 0]*1e3,...
			     'Mp', [167.5 34.6 10.8 0]*1e3,...
			     'M', [185.015 38 12 3.38]*1e3,...
			     'tb',[150 125 600 inf],...
			     'Sref',pi*(5)^2/4,...
                 'guidingTime' , startTimeGuiding,...
			     'g0',9.81);
    steering_path = intMultiStage2DWithDrag([startTimeGT endTimeGT], vertical_path.y(:,end) + [0; 3.08*pi/180; 0; 0], params);%steering angle injected
    params = paramSet('T', [2718 798 83 0]*1e3,...
			     'Mp', [167.5 34.6 10.8 0]*1e3,...
			     'M', [185.015 38 12 3.38]*1e3,...
			     'tb',[150 125 600 inf],...
			     'Sref',pi*(5)^2/4,...
                 'beta0' ,pi/2 - steering_path.y(2,end),...
                 'guidingTime' , startTimeGuiding,...
                 'throttleSwitch', true,...
			     'g0',9.81);    
    % período 2*pi*(Rt+max_height)/max_speed = 6534 s.
    guided_path = intMultiStage2DWithDrag([endTimeGT 875+6534], steering_path.y(:,end), params);


speed_vector=guided_path.y(1,:);
height_vector=guided_path.y(3,:);
max_speed=max(speed_vector);
max_height=max(height_vector);
orbital_speed = sqrt(earth_gravitacional_constant/(Rt + max_height)); %m/s

rinyeccion=max_height+params.Rt;
vinyeccion=max_speed;
SMAteorico=params.Rt+1230e3
SMA=earth_gravitacional_constant/2*(1/(earth_gravitacional_constant/rinyeccion-vinyeccion^2/2));
error_SMA=(1-SMAteorico/SMA)*100
h=rinyeccion*vinyeccion*sin(guided_path.y(2,end));
p_=h^2/earth_gravitacional_constant;
eTeorico=0;
e=sqrt(1-p_/SMA);
error_e=e-eTeorico;
Combustible_sobrante=(875-guided_path.x(end))*params.Mp(3)/params.tb(3);
V_ideal=params.Isp(1)*log(params.r(1))*params.g0+params.Isp(2)*log(params.r(2))*params.g0+params.Isp(3)*log(params.r(3))*params.g0
PerdidasTot=V_ideal-vinyeccion


figure(1);
drawCircle; hold on;
plotOrbit(vertical_path,'.-b');
plotOrbit(steering_path,'.-b');
plotOrbit(guided_path,'.-b');