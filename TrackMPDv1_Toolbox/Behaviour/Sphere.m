classdef Sphere 
% Create a Sphere at coordinates center X and center Y 
% with radius

% Zhiyao, S.,Tingting, W., Fumin, X., Ruijie, L. (2008) A simple formula for 
% predicting settling velocity of sediment particles. Water Science and
% Engineering, Vol. 1, No. 1, 37– 43.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: mp_zhiyao2008 June 2018 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        radius;
    end
    methods
        function obj = Sphere(radius)      
            obj.radius = radius;
            if radius*2>0.005 %IJR 04/06/2018
                disp('Warning: the selected microplastic size is higher than 5 mm')
            end
        end
        function Area = getarea(obj)
            Area = 4*pi*obj.radius^2;
        end
        function Volume = getvolume(obj)
            Volume = 4/3*pi*obj.radius^3;
        end
    end
end

