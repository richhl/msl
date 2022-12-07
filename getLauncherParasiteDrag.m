function drag_coefficient = getLauncherParasiteDrag(mach_number, constant_Reynolds, height )
    arguments
        mach_number double
        constant_Reynolds = true
        height (1,1) double = 0
    end
    % Parasite Drag of the Launcher
    % Drag depends on Mach and Height, but
    % you can force Re = 5e7 (by unit length)
    % passing constant_Reynolds = true
    % % Cd = CD Sumatory for each 'similar shape' geometry. 
    % being the geometries of the launcher:
    % 1. Warhead -Cofia : Cd_wave hipersonic || Cd_shape subsonic + Cd_friction
    % 2. Head Cylinder : Cd_friction
    % 3. Cone trunk : Cd_friction
    % 4. Tail Cylinder: Cd_friction + Cd_base
    
    
    %First going to take into  account transonic mach as it implies
    %interpolation and calling twice this function
    %Transonic drag has to be a linear interpolation between sub-hyper mach
    %points drag values
    if mach_number > 0.9 && mach_number < 1.2
        cd_mach_0_9 = getLauncherParasiteDrag(0.9, constant_Reynolds, height);
        cd_mach_1_2 = getLauncherParasiteDrag(1.2, constant_Reynolds, height);
        drag_coefficient = cd_mach_0_9 + ((cd_mach_1_2-cd_mach_0_9)/(1.2-0.9)) * (mach_number-0.9);
        return
    end
    
        
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
    
    
    
    %%Go first for friction drag as we need it for subsonic shape drag calc 
    
    % raised https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p15316
    % advised to use warhead, cylinder1+cylinder2, cone trunk, cylinder3+cylinder4
    % Friction Drag for warhead
    reference_length = warhead_length;
    wet_surface = warhead_surface; %surface in contact with the flow
    element_shape = shapeType.cone;
    Cd_friction_wh = getFrictionDragCoefficient(mach_number, reference_length, element_shape, wet_surface, reference_surface, constant_Reynolds, height);
    
    
    % Friction Drag for head_cylinder = cylinder1 + cylinder2 (same shape)
    reference_length = cylinder1_length + cylinder2_length;
    wet_surface = pi * warhead_d * reference_length; %surface in contact with the flow
    element_shape = shapeType.flat;
    Cd_friction_head_cylinder = getFrictionDragCoefficient(mach_number, reference_length, element_shape, wet_surface, reference_surface, constant_Reynolds, height);
    
    % Friction Drag for Cone trunk
    reference_length = cone_trunk_length;
    wet_surface = pi * 1/2*(warhead_d + d_base) * sqrt(reference_length^2 + (1/2*(d_base - warhead_d)^2)); %surface in contact with the flow
    element_shape = shapeType.cone;
    Cd_friction_cone_trunk = getFrictionDragCoefficient(mach_number, reference_length, element_shape, wet_surface, reference_surface, constant_Reynolds, height);
    
    % Friction Drag for tail_cylinder = cylinder3 + cylinder4
    reference_length = cylinder3_length + cylinder4_length;
    wet_surface = pi * warhead_d * reference_length; %surface in contact with the flow
    element_shape = shapeType.flat;
    Cd_friction_tail_cylinder = getFrictionDragCoefficient(mach_number, reference_length, element_shape, wet_surface, reference_surface, constant_Reynolds, height);
    
    %Total friction drag is a weigthed sumatory of each 'similar shape geometry' launcher section 
    Cd_friction = Cd_friction_wh + Cd_friction_head_cylinder + ....
                  Cd_friction_cone_trunk + Cd_friction_tail_cylinder;
    
    
    % Subsonic Cd_shape = coefficient * Cd_friction_launcher? question
    % raised https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p15316
    % advised to use  d_base + total_length to get coefficient value
    %% Wave drag for launcher
    theta_cone = atan(warhead_d/(2*warhead_length))*180/pi;
    warhead_surface = pi*warhead_d/2*sqrt(warhead_length^2 + (warhead_d/2)^2);
    geometry.theta_cone = theta_cone;
    Cd_wave = getWaveDrag(mach_number, shapeType.cone, geometry, slenderness, Cd_friction);
    
    
    %% Base drag for launcher
    gas_outlet_surface = 0.8 * reference_surface;
    base_surface = reference_surface;
    Cd_base = getBaseDrag(mach_number, gas_outlet_surface, base_surface, true);
    
    %% Total parasite drag
    drag_coefficient = Cd_wave + Cd_friction + Cd_base;
    return
end
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
                drag_coefficient = 0.25/mach_number * (base_surface - gas_outlet_surface) / base_surface;
        end
    elseif flowKind == flowVelocityRegime.subsonic
        switch engine_on
            case false
                drag_coefficient = (0.12 + 0.13*mach_number^2) * (base_surface / base_surface);
            case true
                drag_coefficient = (0.12 + 0.13*mach_number^2) * (base_surface - gas_outlet_surface) / base_surface;
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
   
   
function drag_coefficient = getFrictionDragCoefficient(mach_number, reference_length, element_shape_type, wet_surface, reference_surface, constant_Reynolds, height)
    mach_number
    reference_length
    constant_Reynolds
    height
    Cd_friction_incompres = getFrictionDragIncompressibleCf(mach_number, reference_length, constant_Reynolds, height);
    Cd_friction_compres = getFrictionDragCompressibleCf(mach_number, reference_length, Cd_friction_incompres, constant_Reynolds, height);
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


function drag_coefficient = getFrictionDragIncompressibleCf(mach_number, reference_length, constant_Reynolds, height)
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
    reynolds_number = getReynoldsNumber(mach_number,reference_length, constant_Reynolds, height);
    if isTurbulentFlow(reynolds_number)
       drag_coefficient = 0.288 * log10( reynolds_number)^-2.45; %Turbulent incompressible
    else
       drag_coefficient = 0.664 / sqrt( reynolds_number); %Laminar incompressible
    end
end

function drag_coefficient = getFrictionDragCompressibleCf(mach_number, reference_length, Cd_friction_incompressible, constant_Reynolds, height)
    gama = 1.4; %adiabatic air constant
    % 4 cases taken into account - Subsonic Laminar/Turbulent, 
    % Supersonic Laminar/Turbulent
    flowKind = getFlowVelocityRegime(mach_number);
    reynolds_number = getReynoldsNumber(mach_number,reference_length, constant_Reynolds, height);
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
function Re = getReynoldsNumber(mach_number, reference_length, constant_Reynolds, height)
    arguments
        mach_number
        reference_length
        constant_Reynolds = false
        height (1,1) double = 0
    end

    % As advised for this work
    % https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p19599
    if constant_Reynolds
        Re = reference_length * 5e7;
        return
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
    mu0 = 3.62e-7 * 100 *0.4536 / 0.3048; %from lb-sec/ft to kg / sec*m)
    T = temperature / 1.8; %kelvin to rankine required by empyrical formula
    T0 = 518.7; %Rankine
    mu = mu0 * ((T / T0)^1.5) * ((T0 + 198.72) / (T + 198.72));
end
