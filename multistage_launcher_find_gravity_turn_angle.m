close all
clc

%%%%%
%% This script aim is finding an initial angle such as the launcher
%% would reach 1230 Km height with enough speed for a circular orbit.
%% (error tolerance: 5% for major axis and excentricity less than 0.03)
%% 2a) solution VA_Ms_practicas_ejecicio1_22-23_enunciado _V3
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

% Our strategy is shooting in a range from 0.35 to 0.45 degrees, and plot
% speed vs. orbital speed. If we got enough speed to leave the launcher
% at 1230 good, otherwise we will select and angle to leave the launcher 
% where speed and orbital speed are almost the same.

%% Uncomment this block if wanting to analyse first shooting iteration
% final_gamma_vector = zeros(1,20);
% final_height_vector = zeros(1,20);
% final_speed_vector = zeros(1,20);
% orbital_speed_vector = zeros(1,20);
% final_theta_vector = zeros(1,20);
% 
% gammaGT=[0:0.01:1];
% for i=30:50
%     vertical_path = intMultiStage2DWithDrag([0 30],[1e-7; 0; 1e-7; 0], params);
%     steering_path = intMultiStage2DWithDrag([30 300], vertical_path.y(:,end) + [0; (gammaGT(i))*pi/180; 0; 0], params);%steering angle injected
%     params = paramSet('T', [2718 798 83 0]*1e3,...
% 			     'Mp', [167.5 34.6 10.8 0]*1e3,...
% 			     'M', [185.015 38 12 3.38]*1e3,...
% 			     'tb',[150 125 600 inf],...
% 			     'Sref',pi*(5)^2/4,...
%                  'beta0' ,pi/2 - steering_path.y(2,end),...
% 			     'g0',9.81);
%     guided_path = intMultiStage2DWithDrag([300 875], steering_path.y(:,end), params);
% 
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

%% First shooting give as gammaGT between 0.38 - 0.39
%% Going to second shooting.

final_gamma_vector = zeros(1,20);
final_height_vector = zeros(1,20);
final_speed_vector = zeros(1,20);
orbital_speed_vector = zeros(1,20);
final_theta_vector = zeros(1,20);

gammaGT=[0.38:5e-4:0.39];
for i=1:length(gammaGT)
    vertical_path = intMultiStage2DWithDrag([0 30],[1e-7; 0; 1e-7; 0], params);
    steering_path = intMultiStage2DWithDrag([30 300], vertical_path.y(:,end) + [0; (gammaGT(i))*pi/180; 0; 0], params);%steering angle injected
    params = paramSet('T', [2718 798 83 0]*1e3,...
			     'Mp', [167.5 34.6 10.8 0]*1e3,...
			     'M', [185.015 38 12 3.38]*1e3,...
			     'tb',[150 125 600 inf],...
			     'Sref',pi*(5)^2/4,...
                 'beta0' ,pi/2 - steering_path.y(2,end),...
			     'g0',9.81);
    guided_path = intMultiStage2DWithDrag([300 875], steering_path.y(:,end), params);

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
txt = '\leftarrow   Selected GammaGT 0.3832'
text(0.3832,7267,txt,'HorizontalAlignment','left')
plot(final_gamma_vector,orbital_speed_vector,'green-');