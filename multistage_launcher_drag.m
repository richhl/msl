%%Parasite Drag of the Launcher
% Cd = CD Sumatory for each 'similar shape' geometry. 
% being the geometries of the launcher:
% 1. Warhead -Cofia : (Cd_shape + Cd_friction)
% 2. Cylinder : (Cd_shape + Cd_friction)
% 3. Cone trunk : (Cd_shape + Cd_friction)
% 4. Tail Cylinder: (Cd_shape + Cd_friction + Cd_base)

%% Global data
warhead_d = 3.5; %m
reference_surface = pi*warhead_d^2/4; %m^2
warhead_length = 4.6;
cylinder1_length = 4.5; %m
cylinder2_length = 10.7; %m
cone_trunk_length = 5; %m cone trunk min diameter = d_warhead, max = d_base
cylinder3_length = 6.6; %m
cylinder4_length = 21; %m
d_base = 5; %m
total_length = warhead_length + cylinder1_length + cylinder2_length...
                + cylinder3_length + cylinder4_length + d_base;
slenderness = total_length/warhead_d; 



clc

%% Cd_warhead = Cd_shape_warhead + Cd_friction_warhead

% Supersonic Cd_wave_warhead - warhead is replaced by a cone as
% required by exercise text

interval = 0.05;
M_ = 0:interval:15; %mach vector points for plotting
number_of_points = size(M_,2);
transonic_start_point = cast((0.9+interval)/interval + 1,"int16"); %M>0.9<1.2 transonic
supersonic_start_point = cast(1.2/interval + 1,"int16"); %M=1.2 supersonic

% Wave drag for launcher
theta_cone = atan(warhead_d/(2*warhead_length))*180/pi;
warhead_surface = pi*warhead_d/2*sqrt(warhead_length^2 + (warhead_d/2)^2);
Cd_wave_vector = zeros(1,number_of_points);

for i=supersonic_start_point:number_of_points
    Cd_wave_vector(i) = (0.083 + 0.096/(M_(i)*M_(i))) * (theta_cone/10)^1.69;
end

% Friction Drag for warhead
Cd_friction_wh_vector = zeros(1,number_of_points);
reference_length = warhead_length;
wet_surface = warhead_surface; %surface in contact with the flow
element_shape = shapeType.cone;
for i=1:number_of_points
    Cd_friction_wh = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_wh_vector(i) = Cd_friction_wh;
end

% Friction Drag for cylinder1
Cd_friction_cylinder1_vector = zeros(1,number_of_points);
reference_length = cylinder1_length;
wet_surface = pi * warhead_d * cylinder1_length; %surface in contact with the flow
element_shape = shapeType.flat;
for i=1:number_of_points
    Cd_friction_cylinder1 = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_cylinder1_vector(i) = Cd_friction_cylinder1;
end

% Friction Drag for cylinder2
Cd_friction_cylinder2_vector = zeros(1,number_of_points);
reference_length = cylinder2_length;
wet_surface = pi * warhead_d * cylinder2_length; %surface in contact with the flow
element_shape = shapeType.flat;
for i=1:number_of_points
    Cd_friction_cylinder2 = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
    Cd_friction_cylinder2_vector(i) = Cd_friction_cylinder2;
end

%check for cylinder1+cylinder2
% Cd_friction_cylinder_vector = zeros(1,number_of_points);
% reference_length = cylinder1_length + cylinder2_length;
% wet_surface = pi * warhead_d * reference_length; %surface in contact with the flow
% element_shape = shapeType.flat;
% for i=1:number_of_points
%     Cd_friction_cylinder = getFrictionDragCoefficient(M_(i), reference_length, element_shape, wet_surface, reference_surface);
%     Cd_friction_cylinder_vector(i) = Cd_friction_cylinder;
% end
% % test = Cd_friction_cylinder1_vector + Cd_friction_cylinder2_vector


% Subsonic Cd_shape = coefficient * Cd_friction_launcher? question
% raised https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p15316
% first approach would use d_warhead + total_length to get coefficient value
% for i=1:transonic_start_point-1
%     Cd_shape =  (60/slenderness^3 + 0.0025*slenderness) * Cd_friction_vector(i);
%     Cd_vector(i) = Cd_shape; 
% end

%% Cd_cylinder1 = Cd_shape_cylinder1 + Cd_friction_cylinder1
% Friction Drag
l_cilynder1 = 4.5 + 10.7 + 6.6;
cyl1_surface = pi*warhead_d*l_cilynder1;
Cd_friction_cyl1_vector = zeros(1,number_of_points);
% Turbulent flow applies as required by exercise text
% for i=1:number_of_points
%     Re_number = (rho0 * M_(i) * a0 * l_cilynder1) / nu;
%     Cd_friction_incompres = 0.288 * log10( Re_number)^-2.45; %Turbulent incompressible
%     if i<transonic_start_point
%         Cd_friction_compres = Cd_friction_incompres * (1/(1+0.08*M_(i)*M_(i))); %Turbulent subsonic compressible
%     elseif i>=supersonic_start_point
%         Cd_friction_compres = Cd_friction_incompres * 1/(1+(gama-1)/2*M_(i)*M_(i))^0.467; %Turbulent supersonic compressible
%     else continue
%     end
%     Cd_friction_mean_cyl1 = 1.25*Cd_friction_compres
%     Cd_friction_cyl1 = Cd_friction_mean_cyl1 * cyl1_surface/ref_surface_geometry1;
%     Cd_friction_cyl1_vector(i) = Cd_friction_cyl1
% end



   
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
function Re = getReynoldsNumber(mach_number, reference_length)
    rho0 = 1.22557; %sea level air density kg/m³
    nu  = 1.8e-5;  %air viscosity
    a0 = 340; %sea level sound speed
    Re  = (rho0 * mach_number * a0 * reference_length) / nu;
end