% Raising velocity (Khatmullina and Isachenko, 2017)
% Isabel Jalón-Rojas 13/11/20

function ws=khatmullina(beh);

n1=length(beh.ParticleDensity);
n2=length(beh.ParticleDequi);

n=max(n1,n2);

if n==n1
    beh.ParticleDequi=beh.ParticleDequi.*ones(size(beh.ParticleDensity));
else
    beh.ParticleDensity=beh.ParticleDensity.*ones(size(beh.ParticleDequi));
end

densi_rel=abs(beh.ParticleDensity-beh.WaterDensity)/beh.WaterDensity;


%D_ast=((g*densi_rel./beh.WaterViscosity).^(1/3)).*beh.ParticleDequi;


g = 9.81*1000; % Gravitational acceleration (mm/s^2)
D = beh.Diameter*1000; %mm
L = beh.ParticleSize*1000; %mm
water_viscosity=beh.WaterViscosity*10^6; %mm^2/s

ws = pi./2./water_viscosity.*densi_rel*g.*D.*L./(55.238*L+12.691); %mm

ws=ws/1000; %m
end
