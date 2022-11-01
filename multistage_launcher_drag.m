%%Parasite Drag of the Launcher
% Cd = CD Sumatory for each 'similar shape' geometry. 
% being the geometries of the launcher:
% 1. Warhead -Cofia : (Cd_shape + Cd_friction)
% 2. Cylinder : (Cd_shape + Cd_friction)
% 3. Cone trunk : (Cd_shape + Cd_friction)
% 4. Tail Cylinder: (Cd_shape + Cd_friction + Cd_base)

%% Global data
rho0 = 1.22557 %sea level air density kg/m³
nu  = 1.8e-5  %air viscosity
a0 = 340 %sea level sound speed
gama = 1.4 %adiabatic air constant
d_warhead = 3.5;
ref_surface_geometry1 = pi*d_warhead^2/4
%% Cd_warhead = Cd_shape_warhead + Cd_friction_warhead

% Subsonic Cd_shape_warhead = coefficient * Cd_friction_cylinder? question
% raised https://moodle.upm.es/titulaciones/oficiales/mod/forum/discuss.php?d=9396#p15316

% Supersonic Cd_wave_warhead - warhead is replaced by a cone as
% required by exercise text

interval = 0.05;
M_ = 0:interval:15; %mach vector points for plotting
number_of_points = size(M_,2);
transonic_start_point = cast((0.9+interval)/interval + 1,"int16"); %M>0.9<1.2 transonic
supersonic_start_point = cast(1.2/interval + 1,"int16"); %M=1.2 supersonic

% Wave drag for launcher warhead
l_warhead = 4.6;
theta_cone = atan(d_warhead/(2*l_warhead))*180/pi;
warhead_surface = pi*d_warhead/2*sqrt(l_warhead^2 + (d_warhead/2)^2);
Cd_wave_wh_vector = zeros(1,number_of_points);

for i=supersonic_start_point:number_of_points
    Cd_wave_wh_vector(i) = (0.083 + 0.096/(M_(i)*M_(i))) * (theta_cone/10)^1.69;
end

% Friction Drag
Cd_friction_wh_vector = zeros(1,number_of_points);
% Turbulent flow applies as required by exercise text
for i=1:number_of_points
    Re_number = (rho0 * M_(i) * a0 * l_warhead) / nu;
    Cd_friction_incompres = 0.288 * log10( Re_number)^-2.45; %Turbulent incompressible
    if i<transonic_start_point
        Cd_friction_compres = Cd_friction_incompres * (1/(1+0.08*M_(i)*M_(i))); %Turbulent subsonic compressible
    elseif i>=supersonic_start_point
        Cd_friction_compres = Cd_friction_incompres * 1/(1+(gama-1)/2*M_(i)*M_(i))^0.467; %Turbulent supersonic compressible
    else continue
    end
    Cd_friction_mean_wh = 1.28*Cd_friction_compres
    Cd_friction_wh = Cd_friction_mean_wh * warhead_surface/ref_surface_geometry1;
    Cd_friction_wh_vector(i) = Cd_friction_wh
end