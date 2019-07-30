function [ x, y, ts ] = particle_track_ode_grid2D(X, Y, U, V, tt, ...
                                                tspan, cc, options, odesolver )
% PARTICLE_TRACK_ODE_GRID  Generates particle tracks from a set of currents
% defined on a rectangular grid using PDE solver
%
% Usage: [x,y,ts] = particle_track_ode_andreas(X,Y,U,V,tt,tspan,cc,options,odesolver )
%
% This function is largely inspired by Bruce Lipphardt's trajectories code.
% The principal difference is the method of interpolating from the grid to
% an arbitrary space time point is different (he uses interpn, I use
% interp2 plus linear interpolation by hand between previous and next
% time of current time).  Also, this assumes that 
%
% Inputs
% ------
% X = x-coordinate of rectangular grid
% Y = Y-coordinate of rectangular grid
%
% U,V = currents at locations in X,Y.  The third dimension should
% be time.  Units should be consistent with those in X and tt.
%
% tt = times at which the currents are defined.  Must have length equal
% to the third dimension of U.
%
% tspan = time period over which to generate tracks.  See ode45 for more
% details on possible tspan formats (some of which might not work with
% this function).
%
% cc = a two column matrix of initial positions for the particles.
%
% options = options for the ODE solver.  See odeset, ode45 and the other
% solvers for details on format and creating.  Defaults to no special options.
%
% odesolver = ode solver function to use.  Defaults to ode45.
%
% OUTPUTS:
%
% x,y = coordinates of particles.  Each column is a different starting
% position (row of cc) and each row is a time in ts.
%
% ts = times of positions in x and y.  A column vector.
%
% NOTE: This function appears to be hopelessly slow in matlab R13 (6.5)
% unless you play with the default tolerances of ode45.  See ODESET
% function and, in particular, the absolute and relative tolerances for
% more information.  You might also want to play with MaxStep, the
% maximum time step to take.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: particle_track_ode_grid.m 001 2018-02-12 10:55:10Z ef $
%        modified by Isabel Jalon-Rojas for 3D version
%
% Copyright (C) 2014-2018 Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist( 'odesolver', 'var' )
  odesolver = 'ode45';
end

if ~exist( 'options', 'var' )
  options = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do tracking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ntraj = size(cc,1);

y0 = [ cc(:,1); cc(:,2) ];

[ts,A] = feval( odesolver, @ptrack_ode_worker_func, tspan, y0, options, ...
		X, Y, U, V, tt );

[x,y] = deal( A(:,1:ntraj), A(:,ntraj+1:end) );

%%%%%%----------------------Subfunctions--------------------%%%%%%%
function f = ptrack_ode_worker_func( T, y, X, Y, ux, uy, tt )
% This function will be used to actually do tracking.
% Assumes that T is a scalar.

% Break up X and Y coordinates.
ntraj = length(y)/2;
[xx,yy] = deal( y(1:ntraj), y(ntraj+1:end) );

% Find timesteps cooresponding to T
d = tt - T;

% If outside time range, return NaNs
if all(d<0) | all(d>0)
  f = repmat(NaN,size(y));
  return
end

% Otherwise, linearly interpolate.
dlt = max( d( d<0 ) );
dgt = min( d( d>=0 ) );
if dgt == 0
  [ux,uy] = deal( ux( :, :, d == 0 ), uy( :, :, d == 0 ) ); 
else  
  ilt = d == dlt;
  igt = d == dgt;
  dd = dgt - dlt;
  clt = dgt / dd;
  cgt = -dlt / dd;
  ux = clt * ux( :, :, ilt ) + cgt * ux( :, :, igt );
  uy = clt * uy( :, :, ilt ) + cgt * uy( :, :, igt );
end

% Then fit to specific points of interest.
nn = ~isnan(xx) & ~isnan(yy);
[UX,UY] = deal( repmat( NaN, size(xx) ) );

UX(nn) = interp2( X, Y, ux, xx(nn), yy(nn), 'linear' );
UY(nn) = interp2( X, Y, uy, xx(nn), yy(nn), 'linear' );
f = [ UX(:); UY(:) ];
