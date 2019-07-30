function [Ws]=mp_zhiyao2008(relative_density,water_viscosity,shape,biofouling_thickness)

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

%----------------------
% DEFINE CONSTANTS
%----------------------
g = 9.81; % Gravitational acceleration (m/s^2)
A = 32.2; % no units
B = 1.17; % no units
n = 1.75; % no units

diameter = 2*(shape.radius+biofouling_thickness);

dstar = (relative_density*g/water_viscosity^2).^(1/3).*diameter;

Ws = (water_viscosity./diameter).*dstar.^3.*((3*A/4)^(2/n)+(3*B/4.*dstar.^3).^(1/n)).^(-n/2);     
end
