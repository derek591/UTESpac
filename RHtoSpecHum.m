function q = RHtoSpecHum(RH,P,T)
% RHtoSpecHum converts RH to q See: http://www.mech.utah.edu/~pardyjak/efd/pottempcalc.pdf
% q [kg/kg], RH [%], P [kPa], T [K]

% constants
l_v = 2.5e6;    %Latent Heat of Vaporization (J/kg)
R_v = 461.5;    %Gas Constant for water (J-K/kg)
T_0 = 273.15;   %Reference Temperature
es0 = 6.11;     %Reference Vapor Pressure at 273.15 (hPa)

%Find Saturated Vapor Pressure
es = es0.*exp((l_v/R_v).*(1/T_0-1./(T)));
e = RH./100.*es; % vapor pressure
q = 0.622.*(e./(P.*10)); %(kg/kg) see Thermo by Cengel Ch. 14

end