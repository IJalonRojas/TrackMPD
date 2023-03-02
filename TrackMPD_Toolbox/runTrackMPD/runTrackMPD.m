% Run TrackMPD in 2D or 3D mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: runTrackMPD July 2019 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas 
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runTrackMPD(conf_name,varargin)

conf=feval(conf_name);
Mode=conf.Traj.Mode;

if strcmpi(Mode,'2D')
    runTrackMPD_2D(conf_name,varargin)
elseif strcmpi(Mode,'3D')
    runTrackMPD_3D(conf_name,varargin)
end

end