function beh = behaviour(beh,tspan,ReleaseTime)

%BEHAVIOR_NEW Calculate the settling velocity over time from the physical
%characteristics and behaviour of particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: behaviour_new.m 001 2022-04-14 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(beh.InitialWsOption,'value') && (strcmpi(beh.Biofouling,'no') || strcmpi(beh.Biofouling,'inWS'))
    
    beh.Ws=beh.Ws0*ones(size(tspan));
    
elseif strcmpi(beh.InitialWsOption,'formulation') && strcmpi(beh.Biofouling,'no')
  
    if strcmpi(beh.WsForm,'waldschlager') 
        beh.Ws0=waldschlager(beh);
    elseif strcmpi(beh.WsForm,'khatmullina')
        beh.Ws0=khatmullina(beh);
    elseif strcmpi(beh.WsForm,'dellino')
        beh.Ws0=dellino(beh);
    elseif strcmpi(beh.WsForm,'zhiyao')
        beh.Ws0=zhiyao(beh);
    end
    
    beh.Ws=beh.Ws0*ones(size(tspan));    
                                        
elseif strcmpi(beh.InitialWsOption,'formulation') && strcmpi(beh.Biofouling,'CteCalculated')   
    
    if strcmpi(beh.ParticleShape,'sphere')
        radius=beh.ParticleSize/2;
        beh.ParticleDensity = beh.ParticleDensity*(radius^3/(radius+beh.BiofoulingThickness)^3)+beh.FilmDensity*(1-(radius^3/(radius+beh.BiofoulingThickness)^3));
    elseif strcmpi(beh.ParticleShape,'fiber') || strcmpi(beh.ParticleShape,'fibre')
        radius=sqrt((beh.ParticleDequi)^3/beh.ParticleSize);
        beh.ParticleDensity = beh.ParticleDensity*(radius^2/(radius+beh.BiofoulingThickness)^2)+beh.FilmDensity*(1-(radius^2/(radius+beh.BiofoulingThickness)^2));
    else
        error('This shape category is for the moment incompatible with the biofouling option CteCalculated')
    end
    
    if strcmpi(beh.WsForm,'waldschlager') 
        beh.Ws0=waldschlager(beh);
    elseif strcmpi(beh.WsForm,'khatmullina')
        beh.Ws0=khatmullina(beh);
    elseif strcmpi(beh.WsForm,'dellino')
        beh.Ws0=dellino(beh);
    elseif strcmpi(beh.WsForm,'zhiyao')
        beh.Ws0=zhiyao(beh);
    end
    
    beh.Ws=beh.Ws0*ones(size(tspan));  
    
    
elseif strcmpi(beh.Biofouling,'Rate')  
   
        beh.ParticleDensity=beh.ParticleDensity+ beh.BiofoulingRate*(tspan-ReleaseTime); %IJR0424
        
        if strcmpi(beh.WsForm,'waldschlager') 
            beh.Ws=waldschlager(beh);
        elseif strcmpi(beh.WsForm,'khatmullina')
            beh.Ws=khatmullina(beh);
        elseif strcmpi(beh.WsForm,'dellino')
            beh.Ws=dellino(beh);
        elseif strcmpi(beh.WsForm,'zhiyao')
            beh.Ws=zhiyao(beh);
        end  
        
elseif strcmpi(beh.InitialWsOption,'rate') && (strcmpi(beh.Biofouling,'no') || strcmpi(beh.Biofouling,'inWS'))
    beh.Ws=beh.Ws0+beh.WsRate*(tspan-ReleaseTime); %IJR0424
else
    error('Options to define settling velocity and biofouling may be incompatibles')
end
                                        
                                        
% Degradation
if strcmpi(beh.Degradation,'Rate')  % just available for waldschlager formulation for the moment
     
    beh.ParticleDequi=beh.ParticleDequi+beh.ParticleDequi*beh.DegradationRate/100*(tspan-ReleaseTime); %IJR0424
    beh.Ws=waldschlager(beh);
end              
    
end

