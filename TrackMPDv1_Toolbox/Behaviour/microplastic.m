classdef microplastic
    %MICROPLASTIC contains all the MicroPlastic Geoometical,Physical, and Biological 
    % Features.
    %
    % Microplastic Physical Features
    %   PolymerType
    %   PolymerAcronym
    %   PolymerDensity
    %
    % Microplastic  Geometrical Features
    %   Category
    %   Shape
    %
    % MicroPlastic Biological Features
    %   BiofoulingThickness
    %   BiofoulingRate
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
               
        % Microplastic Physical Features
        PolymerType
        PolymerAcronym
        PolymerDensity
        
        % Microplastic  Geometrical Features
        Category
        Shape
        
        % MicroPlastic Biological Features
        BiofoulingThickness
        BiofoulingRate
    end
    
    methods
        
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
            possCategory = {'Sphere','Film','Fibre'};
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
        
        function P = set.BiofoulingThickness(P, newBiofoulingThickness)
            if newBiofoulingThickness>=0
                P.BiofoulingThickness = newBiofoulingThickness;
            else
                error('Invalid BiofoulingThickness');
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

