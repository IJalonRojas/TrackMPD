classdef Cylinder 
% Create a Cylinder at coordinates center X and center Y 
% with radius
    properties
        radius;
        height;
    end
    methods
        function obj = Cylinder(radius,height)
%             obj@Shape(centerX,centerY);
            obj.radius = radius;
            obj.height = height;
            if height>0.005 || radius>0.005 %IJR 04/06/2018
                disp('Warning: the selected microplastic size is higher than 5 mm')
            end
        end
        function Area = getarea(obj)
            Area = 2*pi*obj.radius*(obj.radius+obj.height);
        end
        function Volume = getvolume(obj)
            Volume = pi*obj.radius^2*obj.height;
        end
    end
end

