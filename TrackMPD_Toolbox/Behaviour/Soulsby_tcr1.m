% Soulsby formulation: critical shear stress threshold (motion)
% I.Jalon-Rojas 14/09/21
%
% INPUTS:
% ParticleDiameter (mm????)
% ParticleDensity (g/cm3)
%
% OUTPUTS:
% critical shear stress threshold (motion)

function tcr1=Soulsby_tcr1(ParticleDiameter,ParticleDensity,WaterDensity)


g=9.81;
water_viscosity = 10^(-6); % m^2/s
relative_density=abs(ParticleDensity-WaterDensity)/WaterDensity;

Dstar = (relative_density*g/water_viscosity^2).^(1/3).*ParticleDiameter;


Ocr1=0.3/(1+1.2*Dstar)+0.055*(1-exp(-0.02*Dstar));

tcr1=Ocr1*WaterDensity*1000*(relative_density)*g*ParticleDiameter;

end