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
        p = p11*exp(-g0/(R*T11)*(height-h11));
    end
    rho = p/(R*T);
    a = sqrt(1.41*p/rho);
    isa_data = [a,rho, T];
end