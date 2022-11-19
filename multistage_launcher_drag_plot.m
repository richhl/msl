%%Parasite Drag of the Launcher for a set of mach points


%% Mach range to analyse
interval = 0.05;
M_ = 0:interval:15; %mach vector points for plotting
number_of_points = size(M_,2);
Cd = zeros(1,number_of_points)

for i=1:number_of_points
    Cd(i) = getLauncherParasiteDrag(M_(i));
end


cd_figure = plot(M_(10:end), Cd(10:end),'b');
hold on;
title("Cd parásita del lanzador en función de Nº de Mach");
