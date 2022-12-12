%%Parasite Drag of the Launcher for a set of mach points


%% Mach range to analyse
interval = 0.05;
M_ = 0:interval:15; %mach vector points for plotting
number_of_points = size(M_,2);
Cd = zeros(1,number_of_points)

for i=1:number_of_points
    Cd(i) = getLauncherParasiteDrag(M_(i));
end


cd_figure = plot(M_(11:end), Cd(11:end),'b');
set(gca,'color', [0.8 0.8 0.8]);
hold on;
title("Cd parásita del lanzador en función de Nº de Mach");
grid on
xlabel('Número de Mach') 
ylabel('Cd') 
xticks([0.5 0.9 1.2 5 10 15])
xtickangle(45)
yticks([0.3 0.6])
txt = '\leftarrow   Interpolación transónico'
text(1.1,0.3,txt,'HorizontalAlignment','left')
plot(M_(11:end), Cd(11:end),'b');
plot(p.Mach, p.Cd0,'r-');
