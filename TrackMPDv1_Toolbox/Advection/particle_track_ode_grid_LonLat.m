function [ ln, lt, h, ts ] = particle_track_ode_grid_LonLat(Lon, Lat, ZZ, ux, uy, uz, VerticalLayer, ...
                                                  varargin )
% PARTICLE_TRACK_ODE_GRID_LONLAT   Generates particle tracks directly in a
% LonLat coordinate system
%
% Usage: [PLon,PLat,TimeStamps] = particle_track_ode_grid_LonLat(Lon,Lat,U,V,tt, ...
%                                            tspan,LL,options,odesolver )
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
% Inputs
% ------
% Lon,Lat = coordinates of data points.  Data points must lie on a
% regular rectangular Lon,Lat grid.
%
% U,V = currents on Lon,Lat grid.  The third dimension of these matrices
% should be time.  Units MUST be CM/S!
%
% tt = times at which the currents are defined.  Must have length equal
% to the third dimension of U.  Should be in datenum units.
%
% tspan = time period over which to generate tracks.  See ode45 for more
% details on possible formats for tspan. Again, should be in datenum units.
%
% LL = a two column matrix of initial positions for the particles.
% Again, coordinates must be in Lon,Lat.
%
% options = options for the ODE solver.  See odeset, ode45 and the other
% solvers for details on format and creating. Defaults
%
% odesolver = ode solver function to use.  Defaults to ode45.
%
% Outputs
% -------
% PLon,PLat = coordinates of particles.  Each column is a different starting
% position (row of cc) and each row is a time in ts.
%
% ts = times of positions in PLon and PLat.  A column vector.
%
% NOTE: This function appears to be hopelessly slow in matlab R13 (6.5)
% unless you play with the default tolerances of ode45.  See ODESET
% function and, in particular, the absolute and relative tolerances for
% more information.  You might also want to play with MaxStep, the
% maximum time step to take.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	$Id: particle_track_ode_grid_LonLat.m 001 2018-02-12 10:55:10Z ef $
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
% Get XXX/CM - remember angles in True convention
[LonPerCM,LatPerCM] = LonLatPerCM( Lon, Lat );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert values in ux and uy to dLon/dDay and dLat/dDay,
% respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ux = ux .* 60 .* 60 .* 24 .* repmat(LonPerCM,[1,1,size(ux,3)]);
uy = uy .* 60 .* 60 .* 24 .* repmat(LatPerCM,[1,1,size(ux,3)]);
uz = uz ./100 .* 60 .* 60 .* 24; %to dM/dDay

clear LonPerCM LatPerCM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now pass on to particle_track_od_grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(VerticalLayer,'sigma2DepthVar')
    [ln,lt,h,ts] = particle_track_ode_grid_sigmaVar( Lon, Lat, ZZ, ux, uy, uz, varargin{:} );
elseif strcmpi(VerticalLayer,'sigma2DepthCte') || strcmpi(VerticalLayer,'Rectangular')
    [ln,lt,h,ts] = particle_track_ode_grid( Lon, Lat, ZZ, ux, uy, uz, varargin{:} ); 
end


%% Here are some default options for ode that appear to work on large
%% pieces of the earth.
%
%   abs_tol = 1.0e-3; % Not sure about this
%   rel_tol = 1.0e-3; % Not sure about this
%   maxstep = 1/24/2; % 1/2 hour
  
%   options = odeset('RelTol',rel_tol,'AbsTol',abs_tol,'MaxStep',maxstep);
