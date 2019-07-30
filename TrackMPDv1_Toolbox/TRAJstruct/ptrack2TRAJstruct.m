function T = ptrack2TRAJstruct( pfunc, varargin );
% PTRACK2TRAJSTRUCT  Puts output of a particle_track_ode-like function
% into a TRAJ structure.
%
% Usage: TRAJ = ptrack2TRAJstruct( ptrack_function_name, ... )
%
% Inputs
% ------
% ptrack_function_name = name of function to do tracking.  For example,
%                        'particle_track_ode_tri_LonLat'.  It is assumed
%                        that the function returns 3 arguments: Lon, Lat,
%                        TimeStamps (does not really have to be Lon and
%                        Lat for this function to work).
% ... = the rest of the arguments for that tracking function
%
% Outputs
% -------
% TRAJ = a trajectory structure with the resultant data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	$Id: ptrack2TRAJstruct.m 001 2018-02-12 10:55:10Z ef $
%
% Copyright (C) 2014-2018 Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We got the funk
[ln,lt,h,ts] = feval( pfunc, varargin{:} );
[ln,lt,h] = deal(ln',lt',h');

% Awww
T = TRAJstruct( size(ln) );

T.Lon = ln;
T.Lat = lt;
T.Depth = h; %IJR 18/05/18

T.TimeStamp = ts(:)';

% initial positions
T.InitialLonLatDepth = [T.Lon(:,1), T.Lat(:,1) T.Depth(:,1)];

% Calculate final positions
I = repmat( 1:size(T.Lon,2), [size(T.Lon,1),1] );
I( ~isfinite(T.Lon) ) = 0;
I = max( I, [], 2 );

%II = sub2ind( size(T.Lon), (1:size(T.Lon,1))', I ); %IJR 24/05/2018 Error
%when one of the values of I is 0 (when only NaNs in a chunk)

T.FinalLonLatDepth(:,1) = T.Lon(:,end);
T.FinalLonLatDepth(:,2) = T.Lat(:,end);
T.FinalLonLatDepth(:,3) = T.Depth(:,end);

% Duration of tracks
%T.TrajectoryDuration = T.TimeStamp(I)' - T.TimeStamp(1); %IJR 24/05/2018 Error
%when one of the values of I is 0 (when only NaNs in a chunk)
%We can include a loop to solve the problem but not interesting for chunks
%(maybe do a script to join data and calculate the duration there)

% Some metadata
T.ProcessingSteps{end+1} = mfilename;
T.OtherMetadata.(mfilename).ptrack_function_name = pfunc;

