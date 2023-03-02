% Soulsby formulation: critical shear stress threshold (suspension)
% I.Jalón-Rojas 14/09/21
%
% INPUTS:
% ParticleDiameter (mm????)
% ParticleDensity (g/cm3)
%
% OUTPUTS:
% critical shear stress threshold (suspension)

function tcr2=Soulsby_tcr2(ParticleDiameter,ParticleDensity,WaterDensity)


g=9.81;
water_viscosity = 10^(-6); % m^2/s
relative_density=abs(ParticleDensity-WaterDensity)/WaterDensity;

Dstar = (relative_density*g/water_viscosity^2).^(1/3).*ParticleDiameter;


Ocr2=0.3/(1+Dstar)+0.1*(1-exp(-0.05*Dstar));
tcr2=Ocr2*WaterDensity*1000*(relative_density)*g*ParticleDiameter;
end