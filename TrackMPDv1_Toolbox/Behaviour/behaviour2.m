classdef behaviour2
    
    % BEHAVIOUR2 contains all the MPD Geoometical,Physical, and/orBiological 
    % Features needed tro track microplastics with the behaviour 2
    %
    % beh=2: macroplastic (or microplastic with a known ws), plastic density>1 g cm^-3, no biofouling, no degradation
    %
    % Microplastic Physical Features
    %   Ws (Settling Velocity)
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
        BehaviourType = 2
        
        % Microplastic Physical Features
        Ws

        % MicroPlastic Biological Features
        BiofoulingThickness=0;
        BiofoulingRate=0;
        FilmDensity=0;
    end
        
    methods
        function P = set.Ws(P, newWs)
            if newWs>0
                P.Ws = newWs;
            else
              error('Invalid Ws for Behaviour 2');
            end
        end
        
    end
          
end    
   
