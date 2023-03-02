% Raising velocity (Dellino, 2008)
% Isabel Jalón-Rojas 13/11/20

function ws=dellino(beh);

n1=length(beh.ParticleDensity);
n2=length(beh.ParticleDequi);

n=max(n1,n2);

if n==n1
    beh.ParticleDequi=beh.ParticleDequi.*ones(size(beh.ParticleDensity));
else
    beh.ParticleDensity=beh.ParticleDensity.*ones(size(beh.ParticleDequi));
end

rho_f=beh.WaterDensity.*1000; %kg./m3
mu_f=beh.WaterViscosity*1000; %Pa.s

densi_rel=abs(beh.ParticleDensity-beh.WaterDensity)/beh.WaterDensity;
g=9.81; %m2/s
% 
% D_ast=((g*densi_rel./beh.WaterViscosity).^(1/3)).*beh.ParticleDequi;

for i=1:n

%raising
w2=-9999;
w=0.00001;
tol=0.000001;

while abs(w-w2)>tol
    
    Re=w.*beh.ParticleDequi(i).*rho_f./mu_f;

    CD=0.9297./(beh.Phi.^1.6.*Re.^0.0799);
    

    w2=sqrt(4/3*beh.ParticleDequi(i)./CD*densi_rel.*g);
    w=w+0.000001;
    
    if beh.ParticleDensity(i)<beh.WaterDensity
        ws(i)=-w2; 
    else
        ws(i)=w2; 
    end

end

end
