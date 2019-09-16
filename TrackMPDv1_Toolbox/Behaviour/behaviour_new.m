function beh = behaviour_new(beh,tspan)

%BEHAVIOR_NEW Caluclate the settling velocity over time from the physical
%characteristics and behaviour of particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: behaviour_new.m 001 2018-06-18 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Behaviour
%     1: macro or microplastic, plastic density<1 g cm^-3, no biofouling, no degradation
%     2: macroplastic (or microplastic with a known ws), plastic density>1 g cm^-3, no biofouling, no degradation
%     3: microplastic (shape,size,density), plastic density>1 g cm^-3, no biofouling, no degradation
%     4: microplastic (shape,size,density), Biofouling (stationary), no degradation
%     5: microplastic (shape,size,density), Biofouling (non-stationary), no degradation
%     6: microplastic (shape,size,density), No Biofouling, degradation


if beh.BehaviourType==1
    
    beh.Ws=beh.Ws*ones(size(tspan));
    
elseif beh.BehaviourType==2

    
    beh.Ws=beh.Ws*ones(size(tspan)); %(m/s)
    
    
elseif beh.BehaviourType==3
        
    beh.RelativeDensity=(beh.PolymerDensity-beh.WaterDensity)/beh.WaterDensity;

    if  beh.PolymerDensity<=beh.WaterDensity
        Ws=0;
    elseif strcmpi(beh.Category,'sphere') || strcmpi(beh.Category,'shortcylinder')
        Ws = mp_zhiyao2008(beh.RelativeDensity,beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness); 
    elseif strcmpi(beh.Category,'longcylinder')
        Ws = mp_khatmullina2017(beh.RelativeDensity,beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness);
    end
    
    beh.Ws=Ws*ones(size(tspan)); %(m/s)
    
    
elseif beh.BehaviourType==4
    
    if strcmpi(beh.Category,'sphere')
        beh.ParticleDensity = beh.PolymerDensity*(beh.Shape.radius^3/(beh.Shape.radius+beh.BiofoulingThickness)^3)+beh.FilmDensity*(1-(beh.Shape.radius^3/(beh.Shape.radius+beh.BiofoulingThickness)^3));
    elseif strcmpi(beh.Category,'shortcylinder') || strcmpi(beh.Category,'longcylinder')
        beh.ParticleDensity = beh.PolymerDensity*(beh.Shape.radius^2/(beh.Shape.radius+beh.BiofoulingThickness)^2)+beh.FilmDensity*(1-(beh.Shape.radius^2/(beh.Shape.radius+beh.BiofoulingThickness)^2));
    end
         
    beh.RelativeDensity=(beh.ParticleDensity-beh.WaterDensity)/beh.WaterDensity;

    if strcmpi(beh.Category,'sphere') || strcmpi(beh.Category,'shortcylinder')
        Ws = mp_zhiyao2008(beh.RelativeDensity,beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness); 
    elseif strcmpi(beh.Category,'longcylinder')
        Ws = mp_khatmullina2017(beh.RelativeDensity,beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness);
    end
    
    beh.Ws=Ws*ones(size(tspan)); %(m/s)
  
    
elseif beh.BehaviourType==5
    
    beh.BiofoulingThickness=beh.BiofoulingThickness0+ beh.BiofoulingRate*(tspan-tspan(1));
    
    
    if strcmpi(beh.Category,'sphere')
        beh.ParticleDensity = beh.PolymerDensity*(beh.Shape.radius^3./(beh.Shape.radius+beh.BiofoulingThickness).^3)+beh.FilmDensity*(1-(beh.Shape.radius^3./(beh.Shape.radius+beh.BiofoulingThickness).^3));
    elseif strcmpi(beh.Category,'shortcylinder') || strcmpi(beh.Category,'longcylinder')
        beh.ParticleDensity = beh.PolymerDensity*(beh.Shape.radius^2./(beh.Shape.radius+beh.BiofoulingThickness).^2)+beh.FilmDensity*(1-(beh.Shape.radius^2./(beh.Shape.radius+beh.BiofoulingThickness).^2));
    end
         
    for ii=1:length(beh.ParticleDensity)
    
        if beh.ParticleDensity(ii)<=beh.WaterDensity
        
            beh.Ws(ii)=0; %(m/s)     
        else
    
            beh.RelativeDensity(ii)=(beh.ParticleDensity(ii)-beh.WaterDensity)/beh.WaterDensity;

             if strcmpi(beh.Category,'sphere') || strcmpi(beh.Category,'shortcylinder')
                beh.Ws(ii) = mp_zhiyao2008(beh.RelativeDensity(ii),beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness(ii)); 
             elseif strcmpi(beh.Category,'longcylinder')
                beh.Ws(ii) = mp_khatmullina2017(beh.RelativeDensity(ii),beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness(ii));
             end
        end
    end
  
    
elseif beh.BehaviourType==6
     
    
    beh.RelativeDensity=(beh.PolymerDensity-beh.WaterDensity)/beh.WaterDensity;

    if  beh.PolymerDensity<=beh.WaterDensity
        beh.Ws=0;
    elseif strcmpi(beh.Category,'sphere') || strcmpi(beh.Category,'shortcylinder')
        beh.Shape=beh.Shape0;
        beh.Shape.radius=beh.Shape0.radius-beh.Shape0.radius*beh.DegradationRate/100*(tspan-tspan(1));
        beh.Shape.radius(beh.Shape.radius<0)=0;
        Ws = mp_zhiyao2008(beh.RelativeDensity,beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness); 
        Ws(isnan(Ws))=0;
        beh.Ws=Ws;
        
    elseif strcmpi(beh.Category,'longcylinder')
        beh.Shape.radius=beh.Shape0.radius-beh.Shape0.radius*beh.DegradationRate/100*(tspan-tspan(1));
        beh.Shape.height=beh.Shape0.height-beh.Shape0.height*beh.DegradationRate/100*(tspan-tspan(1));
        beh.Shape.radius(beh.Shape.radius<0)=0;
        beh.Shape.height(beh.Shape.height<0)=0;
        Ws = mp_khatmullina2017(beh.RelativeDensity,beh.WaterViscosity,beh.Shape,beh.BiofoulingThickness);
        Ws(isnan(Ws))=0;
        beh.Ws=Ws;
    end

end
    
    
end

