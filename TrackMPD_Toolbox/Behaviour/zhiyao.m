% Raising velocity (Zhiyao, 2008)
% Isabel Jalón-Rojas 13/11/20

function ws=zhiyao(beh);

n1=length(beh.ParticleDensity);
n2=length(beh.ParticleDequi);

n=max(n1,n2);

if n==n1
    beh.ParticleDequi=beh.ParticleDequi.*ones(size(beh.ParticleDensity));
else
    beh.ParticleDensity=beh.ParticleDensity.*ones(size(beh.ParticleDequi));
end

densi_rel=abs(beh.ParticleDensity-beh.WaterDensity)/beh.WaterDensity;
g=9.81; %m2/s

D_ast=((g*densi_rel./beh.WaterViscosity.^2).^(1/3)).*beh.ParticleDequi;

g = 9.81; % Gravitational acceleration (m/s^2)
A = 32.2; % no units
B = 1.17; % no units
n = 1.75; % no units

ws = (beh.WaterViscosity./beh.ParticleDequi).*D_ast.^3.*((3*A/4)^(2/n)+(3*B/4.*D_ast.^3).^(1/n)).^(-n/2);     


end
