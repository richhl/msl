%%Parasite Drag of the Launcher
% Cd = CD Sumatory for each 'similar shape' geometry. 
% being the geometries of the launcher:
% 1. Warhead -Cofia : Cd_wave hipersonic || Cd_shape subsonic + Cd_friction
% 2. Head Cylinder : Cd_friction
% 3. Cone trunk : Cd_friction
% 4. Tail Cylinder: Cd_friction + Cd_base

%% Global data
warhead_d = 3.5; %m
warhead_surface = pi*warhead_d^2/4; %m^2
d_base = 5; %m
reference_surface = pi*d_base^2/4; %m^2
warhead_length = 4.6;
cylinder1_length = 4.5; %m
cylinder2_length = 10.7; %m
cone_trunk_length = 5; %m cone trunk min diameter = d_warhead, max = d_base
cylinder3_length = 6.6; %m
cylinder4_length = 21; %m
total_length = warhead_length + cylinder1_length + cylinder2_length...
                + cylinder3_length + cylinder4_length;
slenderness = total_length/d_base; 

%% Mach range to analyse
interval = 0.05;
M_ = 0:interval:15; %mach vector points for plotting
number_of_points = size(M_,2);
transonic_start_point = cast((0.9+interval)/interval + 1,"int16"); %M>0.9<1.2 transonic
supersonic_start_point = cast(1.2/interval + 1,"int16"); %M=1.2 supersonic

%%Go first for friction drag as we need it for subsonic shape drag calc 

% raised https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p15316
% advised to use warhead, cylinder1+cylinder2, cone trunk, cylinder3+cylinder4
% Friction Drag for warhead
Cd_friction_wh_vector = zeros(1,number_of_points);
reference_length = warhead_length;
wet_surface = warhead_surface; %surface in contact with the flow
element_shape = shapeType.cone;
for i=1:number_of_points
    Cd_friction_wh = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_wh_vector(i) = Cd_friction_wh;
end

% Friction Drag for head_cylinder = cylinder1 + cylinder2 (same shape)
Cd_friction_head_cylinder_vector = zeros(1,number_of_points);
reference_length = cylinder1_length + cylinder2_length;
wet_surface = pi * warhead_d * reference_length; %surface in contact with the flow
element_shape = shapeType.flat;
for i=1:number_of_points
    Cd_friction_head_cylinder = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_head_cylinder_vector(i) = Cd_friction_head_cylinder;
end

% Friction Drag for Cone trunk
Cd_friction_cone_trunk_vector = zeros(1,number_of_points);
reference_length = cone_trunk_length;
wet_surface = pi * 1/2*(warhead_d + d_base) * sqrt(reference_length^2 + (1/2*(d_base - warhead_d)^2)); %surface in contact with the flow
element_shape = shapeType.cone;
for i=1:number_of_points
    Cd_friction_cone_trunk = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_cone_trunk_vector(i) = Cd_friction_cone_trunk;
end

% Friction Drag for tail_cylinder = cylinder3 + cylinder4
Cd_friction_tail_cylinder_vector = zeros(1,number_of_points);
reference_length = cylinder3_length + cylinder4_length;
wet_surface = pi * warhead_d * reference_length; %surface in contact with the flow
element_shape = shapeType.flat;
for i=1:number_of_points
    Cd_friction_tail_cylinder = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_tail_cylinder_vector(i) = Cd_friction_tail_cylinder;
end

%Total friction drag is a weigthed sumatory of each 'similar shape geometry' launcher section 
Cd_friction_vector = zeros(1,number_of_points);
Cd_friction_vector = Cd_friction_wh_vector + Cd_friction_head_cylinder_vector + ....
              Cd_friction_cone_trunk_vector + Cd_friction_tail_cylinder_vector;


% Subsonic Cd_shape = coefficient * Cd_friction_launcher? question
% raised https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p15316
% advised to use  d_base + total_length to get coefficient value
% for i=1:transonic_start_point-1
%     Cd_shape =  (60/slenderness^3 + 0.0025*slenderness) * Cd_friction_vector(i);
%     Cd_vector(i) = Cd_shape; 
% end
%% Wave drag for launcher
theta_cone = atan(warhead_d/(2*warhead_length))*180/pi;
warhead_surface = pi*warhead_d/2*sqrt(warhead_length^2 + (warhead_d/2)^2);
Cd_wave_vector = zeros(1,number_of_points);
for i=1:number_of_points
    geometry.theta_cone = theta_cone
    Cd_wave_vector(i) = getWaveDrag(M_(i), shapeType.cone, geometry, slenderness, Cd_friction_vector(i))
end


%% Base drag for launcher
gas_outlet_surface = 0.8 * reference_surface;
base_surface = reference_surface;
Cd_base_vector = zeros(1,number_of_points);
for i=1:number_of_points
    Cd_base_vector(i) = getBaseDrag(M_(i), gas_outlet_surface, base_surface)
end

%% Total parasite drag
Cd_total = Cd_wave_vector + Cd_friction_vector + Cd_base_vector;
%Cd_total = Cd_wave_vector
%Cd_total = Cd_friction_vector;
%Cd_total = Cd_base_vector;

%Transonic drag has to be a linear interpolation between sub-hyper mach
%points drag values
start_drag = Cd_total(transonic_start_point-1);
end_drag = Cd_total(supersonic_start_point);
number_of_transonic_points = cast(supersonic_start_point - transonic_start_point,'double');
point_number = 1;
for i=1:number_of_transonic_points
    Cd_total(transonic_start_point + (i-1)) = start_drag + ((end_drag-start_drag)/(number_of_transonic_points+1)) * i;
end

%cd_figure = plot(M_(1:19), Cd_total(1:19),'b');
%cd_figure = plot(M_(15:25), Cd_total(15:25),'b');
%cd_figure = plot(M_(20:end), Cd_total(20:end),'b');
cd_figure = plot(M_, Cd_total,'b');
hold on;
title("Cd parásita del lanzador en función de Nº de Mach");



function drag_coefficient= getBaseDrag(mach_number, gas_outlet_surface, base_surface, engine_on)
    arguments
        mach_number
        gas_outlet_surface
        base_surface
        engine_on = false
    end
    flowKind = getFlowVelocityRegime(mach_number);
    if flowVelocityRegime.transonic == flowKind
        drag_coefficient=-1; % 'error value' for transonic
        return 
    end
    
    if flowKind == flowVelocityRegime.supersonic
        switch engine_on
            case false
                drag_coefficient = 0.25/mach_number * (base_surface / base_surface);
            case true
                drag_coefficient = 0.25/mach_number * (base_surface - gas_outlet_surface / base_surface);
        end
    elseif flowKind == flowVelocityRegime.subsonic
        switch engine_on
            case false
                drag_coefficient = (0.12 + 0.13*mach_number^2) * (base_surface / base_surface);
            case true
                drag_coefficient = (0.12 + 0.13*mach_number^2) * (base_surface - gas_outlet_surface / base_surface);
        end
    end
end


function drag_coefficient= getWaveDrag(mach_number, element_shape_type, geometry, slenderness, friction_drag_cf)
% reference_length is mandatory for subsonic/transonic
% friction_drag_cf is mandatory for subsonic/transonic points
     flowKind = getFlowVelocityRegime(mach_number);
    if flowVelocityRegime.transonic == flowKind
        drag_coefficient=-1; % 'error value' for transonic
        return 
    end
    
    if flowKind == flowVelocityRegime.supersonic
        switch element_shape_type
            case shapeType.cone
                drag_coefficient = (0.083 + 0.096/(mach_number^2)) * (geometry.theta_cone/10)^1.69;
            otherwise
                error('NIF. Sorry,but wave drag has been implemented only for cone shape.')
        end
    elseif flowKind == flowVelocityRegime.subsonic
        drag_coefficient =  (60/slenderness^3 + 0.0025*slenderness) * friction_drag_cf;
    end
end
   
   
function drag_coefficient = getFrictionDragCoefficient(mach_number, reference_length, element_shape_type, wet_surface, reference_surface)
    if flowVelocityRegime.transonic == getFlowVelocityRegime(mach_number)
        drag_coefficient=-1; % 'error value' for transonic
        return 
    end
    Cd_friction_incompres = getFrictionDragIncompressibleCf(mach_number, reference_length);
    Cd_friction_compres = getFrictionDragCompressibleCf(mach_number, reference_length, Cd_friction_incompres);
    Cd_friction_mean_wh = getMeanFrictionDrag(element_shape_type, Cd_friction_compres);
    drag_coefficient = Cd_friction_mean_wh * wet_surface/reference_surface;
end

% Turbulent flow applies as required by exercise text. Maybe it would be
% worthwhile to allow activate/deactivate regime velocity check and force regime
% in functions below.

function drag_coefficient = getMeanFrictionDrag(element_shape_type, Cd_friction_compres)
        switch element_shape_type
            case shapeType.flat
                drag_coefficient = 1.25*Cd_friction_compres;
            case shapeType.cone
                drag_coefficient = 1.28*Cd_friction_compres;
        end
end


function drag_coefficient = getFrictionDragIncompressibleCf(mach_number, reference_length)
    % getFrictionDragIncompressibleCf  Returns drag cofficient despite of compressibility effects.
    %   mach_number  Parameter to provide mach number of the element subject to drag.
    %   reference_length Parameter to provide characteristic lentgth of element (f.e aerodinamic chord of an airfoil) 
    %
    %   See also getFrictionDragCompressibleCf.

    %TODO control parameter sanity f.e. not a number passed
    %Watch for singularity
    if mach_number==0
        drag_coefficient = 0;
        return
    end
    reynolds_number = getReynoldsNumber(mach_number,reference_length);
    if isTurbulentFlow(reynolds_number)
       drag_coefficient = 0.288 * log10( reynolds_number)^-2.45; %Turbulent incompressible
    else
       drag_coefficient = 0.664 / sqrt( reynolds_number); %Laminar incompressible
    end
end

function drag_coefficient = getFrictionDragCompressibleCf(mach_number, reference_length, Cd_friction_incompressible)
    gama = 1.4; %adiabatic air constant
    % 4 cases taken into account - Subsonic Laminar/Turbulent, 
    % Supersonic Laminar/Turbulent
    flowKind = getFlowVelocityRegime(mach_number);
    reynolds_number = getReynoldsNumber(mach_number,reference_length);
    if flowKind == flowVelocityRegime.supersonic
        if isTurbulentFlow(reynolds_number) %Turbulent supersonic compressible
            drag_coefficient = Cd_friction_incompressible * 1/(1+(gama-1)/2*mach_number^2)^0.467; 
        else  %Laminar supersonic compressible
            drag_coefficient = Cd_friction_incompressible;
        end
    elseif flowKind == flowVelocityRegime.subsonic
        if isTurbulentFlow(reynolds_number) %Turbulent subsonic compressible
            drag_coefficient = Cd_friction_incompressible * (1/(1+0.08*mach_number^2)); 
        else  %Laminar subsonic compressible
             drag_coefficient = Cd_friction_incompressible * (1/(1+0.17*mach_number^2))^0.1295;            
        end
    end
end

function turbulent = isTurbulentFlow(reynolds_number)
    turbulent = abs(reynolds_number) >= 1e6;
end

function regime = getFlowVelocityRegime(mach_number)
    mach_number = abs(mach_number);
    if mach_number >= 1.2
        regime = flowVelocityRegime.supersonic;
    elseif mach_number <= 0.9
        regime = flowVelocityRegime.subsonic;
    else
        regime = flowVelocityRegime.transonic;
    end
end

%% Would return sea level Reynolds otherwise specified
function Re = getReynoldsNumber(mach_number, reference_length, height)
    arguments
        mach_number
        reference_length
        height (1,1) double = 0
    end
%     rho0 = 1.22557; %sea level air density kg/m³
%     mu  = 1.8e-5;  %air viscosity
%     a0 = 340; %sea level sound speed    isaData = getISAValuesFromHeight(height)
    isaData = getISAValuesFromHeight(height);
    a=isaData(1);
    rho = isaData(2);
    temperature = isaData(3);
    mu = getAirDynamicViscosity(temperature);
    Re  = (rho * mach_number * a * reference_length) / mu;
end

%From sutherland 's formula https://www.grc.nasa.gov/WWW/K-12/airplane/viscosity.html
function mu = getAirDynamicViscosity(temperature)
%Temperature Kelvin
    mu0 = 3.62e-7 * 100 *0.4536 / 0.3048 %from lb-sec/ft to kg / sec*m)
    T = temperature / 1.8 %kelvin to rankine required by empyrical formula
    T0 = 518.7 %Rankine
    mu = mu0 * ((T / T0)^1.5) * ((T0 + 198.72) / (T + 198.72));
end


%inspired by https://github.com/LucianoGP95/International-Standard-Atmosphere-calculator...
% /blob/main/ISA%20calculator/International_Standard_Atmosphere.mlx
function isa_data = getISAValuesFromHeight(height)
    g0= 9.81; %m/s^2  
    To = 288.15; %Kelvin
    po = 101325; %Pa
    R = 287.04; %Air constant
    p11 = 22632; %Pa Above tropopause
    T11 = 216.65; %K Above tropopause
    h11 = 11000; %meters
    a0 = 340.294; %m/s
    
    %Solves air conditions for a specific height
    if height<11001 %Under tropopause
        T = To-6.5*(height/1000);
        p = po*(1-0.0065*(height/To))^5.2561;
    else %Above tropopause
        T = 216.65;
        p = p11*exp(-g/(R*T11)*(height-h11));
    end
    rho = p/(R*T);
    a = sqrt(1.41*p/rho);
    isa_data = [a,rho, T];
end