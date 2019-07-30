function [LonPerM, LatPerM] = LonLatPerM( Lon, Lat, varargin )
% LONLATPERM  Calculates the Longitude and Latitude displacement
% equivalent to 1 m at given locations on the earth.
%
% Usage: [LonPerM, LatPerM] = LonLatPerM( Lon, Lat, spheroid )
%
% This function uses m_map toolbox to make calculations
% 
% Inputs
% ------
% Lon,Lat = equal size matrices of Lon, Lat coordinates
% spheroid = see m_fdist for more details
%
% Outputs
% -------
% LonPerM = change in Lon for 1 cm movement in longitudinal direction at
%            locations in Lon,Lat.
% LatPerM = change in Lat for 1 cm movement in latitudinal direction at
%            locations in Lon,Lat.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: LonLatPerM.m 001 2018-02-12 10:55:10Z ef $
%
% Copyright (C) 2014-2018 Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get XXX/CM - remember angles in True convention
LonPerM = mod( m_fdist( Lon, Lat, 90, 1.0, varargin{:} ) - Lon, 360 );
[tt,ll] = m_fdist( Lon, Lat,  0, 1.0, varargin{:} );
LatPerM = mod( ll - Lat, 360 );
