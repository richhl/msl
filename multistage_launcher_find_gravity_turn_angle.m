close all
clc

%function gammaGT = findGravityTurnAngle(height, start_angle, end_angle)
%%%%%
%% This script aim is finding an initial angle such as the launcher
%% would reach 650 Km height with enough speed for a circular orbit.
%% (error tolerance: 5% for major axis and excentricity less than 0.03)
%% 2a) solution VA_Ms_practicas_ejecicio1_22-23_enunciado_V4
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
endTimeGT = 600; %s time when gravity turn 'ends' 
startTimeGuiding = endTimeGT; % time when guiding by nozzle law enters 
% Our strategy is shooting in a range from 1.5 to 2.5 degrees, and plot
% speed vs. orbital speed. If we got enough speed to leave the launcher
% at 650Km good, otherwise we will select and angle to leave the launcher 
% where speed and orbital speed are almost the same.

%% Uncomment this block if wanting to analyse first shooting iteration
final_gamma_vector = zeros(1,20);
final_height_vector = zeros(1,20);
final_speed_vector = zeros(1,20);
orbital_speed_vector = zeros(1,20);
final_theta_vector = zeros(1,20);

gammaGT=[3.5:-0.05:2];
for i=1:length(gammaGT)
    vertical_path = intMultiStage2DWithDrag([0 startTimeGT],[1e-7; 0; 1e-7; 0], params);
    params = paramSet('T', [2718 798 83 0]*1e3,...
			     'Mp', [167.5 34.6 10.8 0]*1e3,...
			     'M', [185.015 38 12 3.38]*1e3,...
			     'tb',[150 125 600 inf],...
			     'Sref',pi*(5)^2/4,...
                 'guidingTime' , startTimeGuiding,...
			     'g0',9.81);
    steering_path = intMultiStage2DWithDrag([startTimeGT endTimeGT], vertical_path.y(:,end) + [0; (gammaGT(i))*pi/180; 0; 0], params);%steering angle injected
    params = paramSet('T', [2718 798 83 0]*1e3,...
			     'Mp', [167.5 34.6 10.8 0]*1e3,...
			     'M', [185.015 38 12 3.38]*1e3,...
			     'tb',[150 125 600 inf],...
			     'Sref',pi*(5)^2/4,...
                 'beta0' ,pi/2 - steering_path.y(2,end),...
                 'guidingTime' , startTimeGuiding,...
                 'throttleSwitch' , false,...
			     'g0',9.81);    
    guided_path = intMultiStage2DWithDrag([endTimeGT 875], steering_path.y(:,end), params);
    speed_vector=guided_path.y(1,:);
    height_vector=guided_path.y(3,:);
    max_speed=max(speed_vector);
    max_height=max(height_vector);

    final_gamma_vector(i) = gammaGT(i);
    final_beta=90-guided_path.y(2,end)*180/pi;
    final_height_vector(i) = max_height;
    final_speed_vector(i) = max_speed;
    final_beta_vector(i) = final_beta;
    orbital_speed_vector(i) = sqrt(earth_gravitacional_constant/(Rt + max_height)); %m/s
end

%% Plot speed and height vs. orbital speed 
cd_figure = plot(final_gamma_vector,final_speed_vector,'k');
set(gca,'color', [0.8 0.8 0.8]);
hold on;
title("Velocidad vs. velocidad orbital (para distintos ángulos shooting)");
grid on
xlabel('Ángulo de giro por gravedad') 
ylabel('Velocidades') 
xtickangle(45)
plot(final_gamma_vector,orbital_speed_vector,'green-');

% TODO: selectedGamma = interx
%% First shooting give gammaGT 2.02
%% No need to go to second shooting.

% final_gamma_vector = zeros(1,20);
% final_height_vector = zeros(1,20);
% final_speed_vector = zeros(1,20);
% orbital_speed_vector = zeros(1,20);
% final_theta_vector = zeros(1,20);
% 
% gammaGT=[2.015:0.1:2.025];
% for i=1:length(gammaGT)
%     vertical_path = intMultiStage2DWithDrag([0 startTimeGT],[1e-7; 0; 1e-7; 0], params);
%     params = paramSet('T', [2718 798 83 0]*1e3,...
% 			     'Mp', [167.5 34.6 10.8 0]*1e3,...
% 			     'M', [185.015 38 12 3.38]*1e3,...
% 			     'tb',[150 125 600 inf],...
% 			     'Sref',pi*(5)^2/4,...
%                  'guidingTime' , startTimeGuiding,...
% 			     'g0',9.81);
%     steering_path = intMultiStage2DWithDrag([startTimeGT endTimeGT], vertical_path.y(:,end) + [0; (gammaGT(i))*pi/180; 0; 0], params);%steering angle injected
%     params = paramSet('T', [2718 798 83 0]*1e3,...
% 			     'Mp', [167.5 34.6 10.8 0]*1e3,...
% 			     'M', [185.015 38 12 3.38]*1e3,...
% 			     'tb',[150 125 600 inf],...
% 			     'Sref',pi*(5)^2/4,...
%                  'beta0' ,pi/2 - steering_path.y(2,end),...
%                  'guidingTime' , startTimeGuiding,...
% 			     'g0',9.81);    guided_path = intMultiStage2DWithDrag([endTimeGT 875], steering_path.y(:,end), params);
%     speed_vector=guided_path.y(1,:);
%     height_vector=guided_path.y(3,:);
%     max_speed=max(speed_vector);
%     max_height=max(height_vector);
% 
%     final_gamma_vector(i) = gammaGT(i);
%     final_beta=90-guided_path.y(2,end)*180/pi;
%     final_height_vector(i) = max_height;
%     final_speed_vector(i) = max_speed;
%     final_beta_vector(i) = final_beta;
%     orbital_speed_vector(i) = sqrt(earth_gravitacional_constant/(Rt + max_height)); %m/s
% end
% 
% %% Plot speed and height vs. orbital speed 
% cd_figure = plot(final_gamma_vector,final_speed_vector,'k');
% set(gca,'color', [0.8 0.8 0.8]);
% hold on;
% title("Velocidad vs. velocidad orbital (para distintos ángulos shooting)");
% grid on
% xlabel('Ángulo de giro por gravedad') 
% ylabel('Velocidades') 
% xtickangle(45)
% plot(final_gamma_vector,orbital_speed_vector,'green-');
