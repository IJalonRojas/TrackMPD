function [Ws]=mp_khatmullina2017(relative_density,water_viscosity,shape,biofouling_thickness)

%
% Khatmullina and Isachenko. 2017. Settling velocity of microplastic
% particles of regular shapes, Marine Pollution Bulletin, 114, 871-880
%
% All variables in mm (IJR 04/06/2018)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: mp_khatmullina2017 June 2018 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%----------------------
% DEFINE CONSTANTS
%----------------------
g = 9.81*1000; % Gravitational acceleration (mm/s^2)

D = 2*shape.radius*1000+2*biofouling_thickness*1000; %mm
L = shape.height*1000; %mm

water_viscosity=water_viscosity*10^6; %mm^2/s

Ws = pi/2/water_viscosity*relative_density*g.*D.*L./(55.238*L+12.691); %mm
Ws = Ws/1000; %m
end
