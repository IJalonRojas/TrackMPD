function [ x, y, z, ts ] = particle_track_ode_grid(X, Y, ZZ, U, V, W, tt, ...
                                                tspan, cc, options, odesolver)
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

y0 = [ cc(:,1); cc(:,2); cc(:,3) ];

[ts,A] = feval( odesolver, @ptrack_ode_worker_func, tspan, y0, options, ...
		X, Y, ZZ, U, V, W, tt);

[x,y,z] = deal( A(:,1:ntraj), A(:,ntraj+1:2*ntraj), A(:,2*ntraj+1:end));

%%%%%%----------------------Subfunctions--------------------%%%%%%%
function f = ptrack_ode_worker_func( T, y, X, Y, ZZ, ux, uy, uz, tt)
% This function will be used to actually do tracking.
ntraj = length(y)/3;
[xx,yy,zz] = deal( y(1:ntraj), y(ntraj+1:2*ntraj), y(2*ntraj+1:end) );

% Find timesteps corresponding to T in (day)
d = tt - T;

% If outside time range, return NaNs
if all(d<0) || all(d>0)
  return
end

% Otherwise, linearly interpolate.
dlt = max( d( d<0 ) );
dgt = min( d( d>=0 ) );
if dgt == 0
  [ux,uy,uz] = deal( ux( :, :, :, d == 0 ), uy( :, :, :, d == 0 ), uz( :, :, :, d == 0 ));
else 
  ilt = d == dlt;
  igt = d == dgt;
  dd = dgt - dlt;
  clt = dgt / dd;
  cgt = -dlt / dd;
  ux = clt * ux( :, :, :, ilt ) + cgt * ux( :, :, :, igt );
  uy = clt * uy( :, :, :, ilt ) + cgt * uy( :, :, :, igt );
  uz = clt * uz( :, :, :, ilt ) + cgt * uz( :, :, :, igt );
end

% Then fit to specific points of interest (xx,yy,zz).
%
nn = ~isnan(xx) & ~isnan(yy) & ~isnan(zz);
[UX,UY,UZ] = deal( NaN( size(xx) ) );

% find the index of specific points
p=find(nn); np = sum(nn);

% fit the pom velocity field
len_Z=size(ZZ,3); 
uxi=zeros(len_Z,np); uyi=zeros(len_Z,np); uzi=zeros(len_Z,np); hi=zeros(len_Z,np);
 
 
% Interpolate the velocity field at (xx,yy) position for every ZZ
for i=1:len_Z %IJR 21/06/2019
       
    uxi(i,1:np) = interp2( X, Y, ux(:,:,i), xx(nn), yy(nn), 'linear' );
    uyi(i,1:np) = interp2( X, Y, uy(:,:,i), xx(nn), yy(nn), 'linear' );
    uzi(i,1:np) = interp2( X, Y, uz(:,:,i), xx(nn), yy(nn), 'linear' );
    hi(i,1:np) = interp2( X, Y, ZZ(:,:,i), xx(nn), yy(nn), 'linear' );
        
end
    
uxi(isnan(uxi))=0;
uyi(isnan(uyi))=0;
uzi(isnan(uzi))=0;
    
% Interpolate the velocity field at zz vertical position
for itraj = 1:np   
    if isnan(sum(hi(:,itraj))) | hi(len_Z,itraj)>zz(nn(itraj)) %out of domain or sinking or beaching
  
        UX(p(itraj))=0;
        UY(p(itraj))=0;
        UZ(p(itraj))=0;
    
    else    

        UX(p(itraj))=interp1(hi(:,itraj),uxi(:,itraj),zz(nn(itraj)));
        UY(p(itraj))=interp1(hi(:,itraj),uyi(:,itraj),zz(nn(itraj)));
        UZ(p(itraj))=interp1(hi(:,itraj),uzi(:,itraj),zz(nn(itraj)));

    end  
    
end

f = [ UX(:); UY(:); UZ(:) ];

