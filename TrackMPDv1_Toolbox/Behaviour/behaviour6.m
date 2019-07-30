classdef behaviour6
    %% 
    %BEHAVIOUR6 contains all the MPD Geoometical,Physical, and/orBiological 
    % Features needed tro track microplastics with the behaviour 6
    %
    % beh=6:  microplastic (shape,size,density), plastic density>1 g cm^-3, no biofouling, no degradation
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
    %
    % Microplastic  Geometrical Features
    %   Category
    %   Shape0 (Initial shape)
    %   Shape
    %
    % MicroPlastic Biological Features
    %   BiofoulingThickness
    %   BiofoulingRate
    %   FilmDensity
    %
    % Plastic degradation
    %   DegradationRate
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
        BehaviourType = 3
        
        % Water properties
        WaterDensity
        WaterViscosity
        
        % Microplastic Physical Features
        PolymerType
        PolymerAcronym
        PolymerDensity
        Ws
        RelativeDensity
        
        % Microplastic  Geometrical Features
        Category
        Shape0
        Shape
        
        % MicroPlastic Biological Features
        BiofoulingThickness=0;
        BiofoulingRate=0;
        FilmDensity=0;
        
        % Plastic degradation
        DegradationRate
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
        
        function P = set.Shape0(P, newShape0)
            % Control allowable Shape values
            possShape = {'Sphere','Cylinder'};
            switch class( newShape0 )
                case possShape
                    P.Shape0 = newShape0;
                otherwise
                    error('Invalid Shape');
            end
        end
        
        function P = set.Ws(P, newWs)
            if newWs>=0
                P.Ws = newWs;
            else
              error('Invalid Ws for Behaviour 4');
            end
        end
        
        function P = set.RelativeDensity(P, newRelativeDensity)
            if newRelativeDensity>=0
                P.RelativeDensity = newRelativeDensity;
            else
              error('Invalid RelativeDensity');
            end
        end
        
        
    end
    
end

