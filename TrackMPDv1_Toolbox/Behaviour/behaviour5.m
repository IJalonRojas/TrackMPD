classdef behaviour5
    %% 
    %BEHAVIOUR5 contains all the MPD Geoometical,Physical, and/orBiological 
    % Features needed tro track microplastics with the behaviour 5
    %
    % beh=5: microplastic (shape,size,density), Biofouling (stationary), no degradation
    %
    % Water properties
    %   WaterDensity
    %   WaterViscosity
    %
    % Microplastic Physical Features
    %   PolymerType
    %   PolymerAcronym
    %   PolymerDensity
    %   RelativeDensity
    %   Ws (Settling Velocity)
    %   ParticleDensity
    %
    % Microplastic  Geometrical Features
    %   Category
    %   Shape
    %
    % MicroPlastic Biological Features
    %   BiofoulingThickness0
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
        BehaviourType = 5
        
        % Water properties
        WaterDensity
        WaterViscosity
        
        % Microplastic Physical Features
        PolymerType
        PolymerAcronym
        PolymerDensity
        Ws
        RelativeDensity
        ParticleDensity
        
        % Microplastic  Geometrical Features
        Category
        Shape
        
        % MicroPlastic Biological Features
        BiofoulingThickness0
        BiofoulingThickness
        BiofoulingRate
        FilmDensity
    end
    
    methods
        
        function P = set.WaterDensity(P, newWaterDensity)
            if newWaterDensity>=0.9 && newWaterDensity<=1.05
                P.WaterDensity = newWaterDensity;
            else
                error('Invalid Water Density');
            end
        end
        
        function P = set.WaterViscosity(P, newWaterViscosity)
            if newWaterViscosity>=0 && newWaterViscosity<=2*10^(-6)
                P.WaterViscosity = newWaterViscosity;
            else
                error('Invalid Water Viscosity');
            end
        end
        
        function P = set.PolymerType ( P, newPolymerType)
            if ischar(newPolymerType) && ndims(newPolymerType)==2 && size(newPolymerType,1)==1
                P.PolymerType= newPolymerType;
            else
                error('Invalid PolymerType');
            end
        end
        
        function P = set.PolymerAcronym ( P, newPolymerAcronym)
            if ischar(newPolymerAcronym) && ndims(newPolymerAcronym)==2 && size(newPolymerAcronym,1)==1
                P.PolymerAcronym= newPolymerAcronym;
            else
                error('Invalid PolymerAcronym');
            end
        end
        
        function P = set.PolymerDensity(P, newPolymerDensity)
            if newPolymerDensity>=0
                P.PolymerDensity = newPolymerDensity;
            else
                error('Invalid PolymerDensity');
            end
        end
        
        function P = set.Category(P, newCategory)
            % Control allowable Category values
            possCategory = {'Sphere','ShortCylinder','LongCylinder'}; %IJR: To be discussed
            switch lower(newCategory)
                case lower(possCategory)
                    P.Category = newCategory;
                otherwise
                    error('Invalid Category');
            end
        end
        
        function P = set.Shape(P, newShape)
            % Control allowable Shape values
            possShape = {'Sphere','Cylinder'};
            switch class( newShape )
                case possShape
                    P.Shape = newShape;
                otherwise
                    error('Invalid Shape');
            end
        end
        
        function P = set.Ws(P, newWs)
            if newWs>=0
                P.Ws = newWs;
            else
              error('Invalid Ws for Behaviour 5');
            end
        end
        
        function P = set.ParticleDensity(P, newParticleDensity)
            if newParticleDensity>=0
                P.ParticleDensity = newParticleDensity;
            else
              error('Invalid ParticleDensity');
            end
        end
        
        function P = set.RelativeDensity(P, newRelativeDensity)
            if newRelativeDensity>=0
                P.RelativeDensity = newRelativeDensity;
            else
              error('Invalid RelativeDensity');
            end
        end
        
        function P = set.BiofoulingThickness0(P, newBiofoulingThickness0)
            if newBiofoulingThickness0>=0
                P.BiofoulingThickness0 = newBiofoulingThickness0;
            else
                error('Invalid BiofoulingThickness0 for Behaviour 5');
            end
        end
        
        function P = set.BiofoulingThickness(P, newBiofoulingThickness)
            if newBiofoulingThickness>=0
                P.BiofoulingThickness = newBiofoulingThickness;
            else
                error('Invalid BiofoulingThickness0 for Behaviour 5');
            end
        end
        
        
         function P = set.FilmDensity(P, newFilmDensity)
            if newFilmDensity>0
                P.FilmDensity = newFilmDensity;
            else
                error('Invalid Film Density for Behaviour 5');
            end
        end
        
        function P = set.BiofoulingRate(P, newBiofoulingRate)
            if newBiofoulingRate>=0
                P.BiofoulingRate = newBiofoulingRate;
            else
                error('Invalid BiofoulingRate');
            end
        end
        
    end
    
end

