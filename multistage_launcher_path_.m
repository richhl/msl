close all
clc
%% Trabajo 1
%Parametros
P.g0=9.81
P.Rt=6378e3
P.mu=P.g0*P.Rt^2
P.M=[185015 38000 12000 3380]
P.Mp=[167500 34600 10800 0]
P.Isp=[248.1172 293.8785 470.0419 0]
P.tb=[150 125 600 Inf]
P.M0=fliplr(cumsum(fliplr(P.M)));
P.t0=[0,cumsum(P.tb(1:end-1))];
P.S=1-P.Mp./P.M;
P.r=P.M0./(P.M0-P.Mp)
P.i =1:length(P.t0);
P.thrust=[2718 798 83 0]*1000;
P.Sref = pi*5^2/4;
tGT=30;
%Ascenso vertical 0<t<tGT
tinterval = [0 tGT];
ICVV = [1e-7; 0; 1e-7; 0]; % [m/s; rad; m; rad]	

options_VV = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
solVV = ode45(@(t,y) VueloVertical(t, y, P), tinterval, ICVV, options_VV);


%SOLUCIONES
%t1=solVV.x;
vectorV=solVV.y(1,:);
vectorZ=solVV.y(3,:);

zGT=max(vectorZ) %Altura final
vGT=max(vectorV) %velocidad final    

%Giro por gravedad

ICGG = [vGT; 0; zGT; 0]; % [m/2; rad; m; rad]

%

function yp=VueloVertical(t,y,P)
    	I = interp1(P.t0, P.i, t,'previous', 'extrap');
        a0 = 340 %m/s
	    %%%%% Thrust
	    thrust = (t<P.tb(I))*9.8*P.Isp(I)*P.Mp(I)/P.tb(I)./(P.M0(I)-P.Mp(I).*(t-P.t0(I))/P.tb(I));
	    Drag = (I==1 && t>30)*-0.5*P.Sref*getLauncherParasiteDrag(y(1)/a0)*rhoISA(y(3)).*y(1).*abs(y(1))./(P.M0(I)-P.Mp(I).*(t-P.t0(I))/P.tb(I));
	    %%%%%
    
	    yp = [thrust + Drag - P.g0*cos(y(2))./(1+y(3)/P.Rt).^2;...
              (P.g0./y(1)./(1+y(3)/P.Rt).^2-y(1)./(P.Rt+y(3))).*sin(y(2));...
		      y(1).*cos(y(2));...
		      y(1).*sin(y(2))./(P.Rt+y(3))];
end

function rho = rhoISA(z)
	rho = (1.225*(1-z/44.338e3).^4.25).*(z < 11e3) + (0.365*exp((11e3-z)/(6.35e3))).*(z >= 11e3);
end