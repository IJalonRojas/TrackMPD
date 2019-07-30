classdef behaviour1
    
    % BEHAVIOUR1 contains all the MPD Geoometical,Physical, and/orBiological 
    % Features needed tro track microplastics with the behaviour 1
    %
    % beh=1: plastic density<1 g cm^-3, no biofouling, no degradation
    %
    % Microplastic Physical Features
    %   Ws (Settling Velocity)
    %   PolymerDensity (useful for windage)
    %
    % MicroPlastic Biological Features
    %   BiofoulingThickness
    %   BiofoulingRate
    %   FilmDensity
    %   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % 	$Id: behaviour5 001 2018-06-18 Z ijalonrojas $
    %
    % Copyright (C) 2017-2019 Isabel Jalon-Rojas
    % Licence: GPL (Gnu Public License)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    
    properties
        % Basics
        CreationInfo = 'UNSW'
        CreateTimeStamp = now
        CreateTimeZone = 'GMT'
        BehaviourType = 1
        
        % Microplastic Physical Features
        Ws=0;
        PolymerDensity

        % MicroPlastic Biological Features
        BiofoulingThickness=0;
        BiofoulingRate=0;
        FilmDensity=0;
    end
        
    
     methods
        
        
        function P = set.PolymerDensity(P, newPolymerDensity)
            if newPolymerDensity>=0
                P.PolymerDensity = newPolymerDensity;
            else
                error('Invalid PolymerDensity');
            end
        end
     end
    
end    
   
