% Raising velocity (Waldschlager and Schuttrumpf, 2019)
% Isabel Jalón-Rojas 13/11/20

function ws=waldschlager(beh)

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

D_ast=((g*densi_rel./beh.WaterViscosity).^(1/3)).*beh.ParticleDequi;

for i=1:n

%raising
w2=-9999;
w=0.00001;
tol=0.000001;

while abs(w-w2)>tol
    
    Re=w.*beh.ParticleDequi(i)./beh.WaterViscosity;
    
    if beh.ParticleDensity(i)<beh.WaterDensity && (~strcmpi(beh.ParticleShape,'fiber') && ~strcmpi(beh.ParticleShape,'fibre')) 
        CD=(20./Re+10./sqrt(Re)+sqrt(1.195-beh.CSF)).*(6/beh.P).^(1-beh.CSF);
    
    elseif beh.ParticleDensity(i)<beh.WaterDensity && (strcmpi(beh.ParticleShape,'fiber') || strcmpi(beh.ParticleShape,'fibre')) 
        CD=10./sqrt(Re)+sqrt(beh.CSF);
        
    elseif beh.ParticleDensity(i)>beh.WaterDensity && (~strcmpi(beh.ParticleShape,'fiber') && ~strcmpi(beh.ParticleShape,'fibre')) 
        CD=3./(beh.CSF.*Re^(1/3));
    elseif beh.ParticleDensity(i)>beh.WaterDensity && (strcmpi(beh.ParticleShape,'fiber') || strcmpi(beh.ParticleShape,'fibre')) 
        CD=4.7./sqrt(Re)+sqrt(beh.CSF);
    end
    

    w2=sqrt(4/3*beh.ParticleDequi(i)./CD.*densi_rel(i)*g);
    w=w+0.000001;
    
    if beh.ParticleDensity(i)<beh.WaterDensity
        ws(i)=-w2; 
    else
        ws(i)=w2; 
    end
end




end
