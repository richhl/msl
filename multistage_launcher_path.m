close all
clc

addpath('/home/rich/grado/asignaturas/misiles/practicas/sesion_3_lanzadores/');

%% Trabajo 1
%Parametros
params = paramSet('Isp', [248.1172 293.8785 470.0419 0],...
				 'Mp', [167.5 34.6 10.8 0]*1e3,...
				 'M', [185.015 38 12 3.38]*1e3,...
				 'tb',[150 125 600 inf],...
				 'Sref',pi*(5)^2/4,...
				 'CD0',0,...
				 'g0',9.81);

vertical_path = intMultiStage2DWithDrag([0 30],[1e-7; 0; 1e-7; 0], params);
steering_path = intMultiStage2DWithDrag([30 300], vertical_path.y(:,end) + [0; .45*pi/180; 0; 0], params);%steering angle injected
guided_path = intMultiStage2DWithDrag([300 875], steering_path.y(:,end), params);


speed_vector=guided_path.y(1,:);
height_vector=guided_path.y(3,:);

max_speed=max(speed_vector);
max_height=max(height_vector)


figure
drawCircle; hold on;
plotOrbit(vertical_path,'.-b');
plotOrbit(steering_path,'.-b');
plotOrbit(guided_path,'.-b');


% figure
% subplot(2,1,1);
% plot(steering_path.x,steering_path.y(3,:),'.-k'); grid; hold on;
% xline(steering_path.p.t0(2), 'r')
% xline(steering_path.p.t0(3), 'r')
% xline(steering_path.p.t0(4), 'r')
% %xlim([0 25000])
% ylabel('z (km)')
% xlabel('t (s)')
% subplot(2,1,2);
% plot(steering_path.x,steering_path.y(1,:)./1000,'.-k'); grid; hold on;
% xline(steering_path.p.t0(2), 'r')
% xline(steering_path.p.t0(3), 'r')
% xline(steering_path.p.t0(4), 'r')
% %xlim([0 25000])
% ylabel('v (km/s)')
% xlabel('t (s)')

% figure
% subplot(2,1,1);
% plot(vertical_path.x,vertical_path.y(3,:),'.-k'); grid; hold on;
% xline(vertical_path.p.t0(2), 'r')
% xline(vertical_path.p.t0(3), 'r')
% xline(vertical_path.p.t0(4), 'r')
% xlim([0 25000])
% ylabel('z (km)')
% xlabel('t (s)')
% subplot(2,1,2);
% plot(vertical_path.x,vertical_path.y(1,:)./1000,'.-k'); grid; hold on;
% xline(vertical_path.p.t0(2), 'r')
% xline(vertical_path.p.t0(3), 'r')
% xline(vertical_path.p.t0(4), 'r')
% xlim([0 25000])
% ylabel('v (km/s)')
% xlabel('t (s)')


% P.g0=9.81
% P.Rt=6378e3
% P.mu=P.g0*P.Rt^2
% P.M=[185015 38000 12000 3380]
% P.Mp=[167500 34600 10800 0]
% P.Isp=[248.1172 293.8785 470.0419 0]
% P.tb=[150 125 600 Inf]
% P.M0=fliplr(cumsum(fliplr(P.M)));
% P.t0=[0,cumsum(P.tb(1:end-1))];
% P.S=1-P.Mp./P.M;
% P.r=P.M0./(P.M0-P.Mp)
% P.i =1:length(P.t0);
% P.thrust=[2718 798 83 0]*1000;
% P.Sref = pi*5^2/4;
% 
% tGT=30;
% tSteering=300;
% 
% %Ascenso vertical 0<t<tGT
% tinterval = [0 150];
% ICVV = [1e-7; 0; 1e-7; 0]; % [m/s; rad; m; rad]	
% 
% options_VV = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
% solVV = ode45(@(t,y) Ecs2D(t, y, P), tinterval, ICVV, options_VV);
% 
% %SOLUCIONES
% vectorV=solVV.y(1,:);
% vectorZ=solVV.y(3,:);
% 
% zGT=max(vectorZ) %Altura final
% vGT=max(vectorV) %velocidad final    
% 
% %Giro por gravedad
% % tintervalGG= [tGT tSteering];
% % gammaGT=0*pi/180;
% % ICGG = [vGT; gammaGT; zGT; 0];% [m/2; rad; m; rad]
% % 
% % options_GG = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
% % solGG = ode45(@(t,y) Ecs2D(t, y, P), tintervalGG, ICGG, options_GG);
% % 
% % vectorV_GG=solGG.y(1,:);
% % vectorgamma_GG=solGG.y(2,:);
% % vectorZ_GG=solGG.y(3,:);
% % 
% % zGGfinal=max(vectorZ_GG)
% % gammaGGfinal=max(vectorgamma_GG)
% % vGGfinal=max(vectorV_GG)
% 
%  
% %Maniobra guiado
% 
% %ICMG = [; ; ;]
% 
% function yp=Ecs2D(t,y,P)
%     	I = interp1(P.t0, P.i, t,'previous', 'extrap');
%         ISAValues=getISAValuesFromHeight(y(3));
% 	    a=ISAValues(1);
%         rho=ISAValues(2);
%         T=ISAValues(3);
%         Mach=y(1)/a;
%         %%%%% Thrust
%         %P.thrust(I)
%         %thrust = P.thrust(I)./(P.M0(I)-P.Mp(I).*(t-P.t0(I))/P.tb(I))
%         thrust = (P.g0*P.Isp(I)*P.e(I)/P.tb(I))./(1-P.e(I)*t/P.tb(I))
% 	if (false) %Despreciamos resistencia para Mach<0.5 y para etapas distintas a la 1
%     	%if (I==1 && Mach>=0.5) %Despreciamos resistencia para Mach<0.5 y para etapas distintas a la 1
%             Drag = -0.5*P.Sref*getLauncherParasiteDrag(Mach)*rho.*y(1).*abs(y(1))./(P.M0(I)-P.Mp(I).*(t-P.t0(I))/P.tb(I))
%         else
%             Drag = 0
%         end
%             
% 	    yp = [thrust + Drag - P.g0*cos(y(2))./(1+y(3)/P.Rt).^2;...
%               (P.g0./y(1)./(1+y(3)/P.Rt).^2-y(1)./(P.Rt+y(3))).*sin(y(2));...
% 		      y(1).*cos(y(2));...
% 		      y(1).*sin(y(2))./(P.Rt+y(3))];
% end