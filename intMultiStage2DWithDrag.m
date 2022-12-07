function m = intMultiStage2DWithDrag(tspan, IC, varargin)
	if nargin < 3
		params = paramSet();
	else
		params = varargin{1};
	end

	options = odeset('Events', @EventsFcn, 'RelTol', 1e-10, 'AbsTol', 1e-10);
	sol = ode45(@(t,y) odefun(t, y, params), tspan, IC, options);
	m = sol;
	m.p = params;
	
%%%%%%%%%%%%%%%%%%%Extra, useful quantities
	m.Pd = .5*rhoISA(m.y(3,:)).*m.y(1,:).^2;
	m.Qwp = .5*rhoISA(m.y(3,:)).*m.y(1,:).^3;
	I = interp1(m.p.t0, m.p.i, m.x, 'previous', 'extrap');
	m.Thrust = 9.8*m.p.Isp(I).*m.p.Mp(I)./params.tb(I);
	m.rho = rhoISA(m.y(3,:));
%%%%%%%%%%%%%%%%%%%
end

function yp = odefun(t, y, P)
	I = interp1(P.t0, P.i, t,'previous', 'extrap');
	
    ISAValues = getISAValuesFromHeight(y(3));
    sound_speed = ISAValues(1);
    current_mach = y(1)/sound_speed;

	%%%%% Thrust
	thrust = 9.8*P.Isp(I)*P.Mp(I)/P.tb(I)./(P.M0(I)-P.Mp(I).*(t-P.t0(I))/P.tb(I));
    %Drag taken into account only for stage 1 and mach>=0.5 as fixed for
%     %this problem.
     if (current_mach >= 0.5 && I==1)
         Drag = -0.5*P.Sref*getLauncherParasiteDrag(current_mach)*rhoISA(y(3)).*y(1).*abs(y(1))./(P.M0(I)-P.Mp(I).*(t-P.t0(I))/P.tb(I));
     else
        Drag = 0;
    end
	%%%%%

    %%Guide law tan(beta) = (1 − t/tf) * tan(beta0).
      
    if (t>300 && not(P.beta0==0)) %P.beta0==0 if this function is called for steering path
        beta0 = P.beta0; %Grab beta0 angle from last path.
        tf = P.tb(3) + P.tb(1) + P.tb(2); %final time of guiding (absolute) 
        error = pi/2 - y(2) - atan(((1 - (t)/(tf))*tan(beta0)));
        delta = error;
    else
        delta = 0; %engine line aligned with launcher axis
    end

	yp = [thrust.*cos(delta) + Drag - P.g0*cos(y(2))./(1+y(3)/P.Rt).^2;...
      thrust.*sin(delta)+(P.g0./y(1)./(1+y(3)/P.Rt).^2-y(1)./(P.Rt+y(3))).*sin(y(2));...
	  y(1).*cos(y(2));...
	  y(1).*sin(y(2))./(P.Rt+y(3))];
end

function rho = rhoISA(z)
	rho = (1.225*(1-z/44.338e3).^4.25).*(z < 11e3) + (0.365*exp((11e3-z)/(6.35e3))).*(z >= 11e3);
end

function [position, isterminal, direction] = EventsFcn(t, y)
	position = [y(1).*cos(y(2)); y(3)];
	isterminal = [0; 1];
	direction = [-1; -1];
end
