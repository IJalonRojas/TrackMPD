function [ux,uy] = transformUScalar(Lon,Lat,ux,uy,Coordinates,Direction)
% transformUScalar: 
% Direction == -1 => Lon-Lat/Day to m/s
% Direction == +1 => m/s => Lon-Lat/Day or cm/s
%
% Usage: [ux,uy] = transformUScalar(Lon,Lat,ux,uy,Coordinates,Direction)
%
% This function is essentially the same as the particle_track_ode_grid
% function, except that it converts currents from cm/s to dLon/dDay and
% dLat/dDay and does particle tracking directly in the Lon,Lat coordinate
% system.  This has the advantage that it should work fairly seamlessly over
% large areas on a spherical earth so long as the individual pieces of the
% track are small enough that the earth can be assumed locally flat and the
% dLon and dLat make sense over those scales.  The size of individual pieces
% of the track can be set sufficiently small for most parts of the earth
% (i.e., all but close to the poles) by choosing appropriate values of the
% options argument.
%
% Also, it is assumed in this function that the separation between nearby
% grid points is small enough that straight LonLat lines between grid points
% are not significantly different from great arcs on the surface of the
% earth.
%
% This function uses the m_idist function to estimate the number of cm
% per degree of longitude and latitude, so the M_MAP toolbox must be on
% the path.
%
% IMPORTANT NOTE: This function has a fundamental problem tracking things
% across the longitudinal boundary (i.e. where it passes from -180 to
% +180, for example).  Therefore, make sure that longitude is continuous
% over the longitudinal range of your data by adding an appropriate
% factor of 360 degrees where necessary.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	$Id: transformUScalar.m 001 2018-02-12 10:55:10Z ef $
%        modified by Isabel Jalon-Rojas for 3D version
%
% Copyright (C) 2014-2018 Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First need to calculate number of km per Lon o Lat
% at locations of data.  Currently uses m_map toolbox
% to do that, but there are probably simple ways to
% do this without the toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(Coordinates,'spherical') 

    %Get XXX/M - remember angles in True convention
    [LonPerCM,LatPerCM] = LonLatPerM( Lon, Lat );

    % Convert ux, uy from m/s to Lon/day and Lat/Day
    if Direction > 0
      ux = ux * 86400 * LonPerCM; 
      uy = uy * 86400 * LatPerCM;
    
    % Convert ux, uy from Lon/day and Lat/Day to m/s 
    else
      ux = ux / 86400 / LonPerCM; 
      uy = uy / 86400 / LonPerCM; 
    end

elseif strcmpi(Coordinates,'cartesian')
    
    % Convert ux, uy from m/s to m/day
    if Direction > 0
      ux = ux*86400; 
      uy = uy*86400;
      
    % Convert ux, uy from m/day to m/s
    else
      ux = ux/86400; 
      uy = uy/86400; 
    end
      
else
    error('The coordinate system should be spherical or cartesian, check the input file')
end


end

