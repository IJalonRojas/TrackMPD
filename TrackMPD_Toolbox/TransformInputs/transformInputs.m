function transformInputs(conf_name,confOGCM_name)

% TRANSFORMOUTPUTS select the function to transmorm the format of OGCM inputs
% to TrackMPD format
%
% INPUTS: model name and mode (2D or 3D)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: transformFVCOMoutputs July 2019 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas 
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf=feval(conf_name);

OGCMmodel_name=conf.OGCM.Model_name;
Mode=conf.Traj.Mode;

if strcmpi(OGCMmodel_name,'POM') && strcmpi(Mode,'3D')
    transformPOMinputs_3D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'POM') && strcmpi(Mode,'2D')
    transformPOMinputs_2D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'FVCOM') && strcmpi(Mode,'3D')
    transformFVCOMinputs_3D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'FVCOM') && strcmpi(Mode,'2D')
    transformFVCOMinputs_2D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'TELEMAC') && strcmpi(Mode,'3D')
    transformTELEMACinputs_3D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'TELEMAC') && strcmpi(Mode,'2D')
    transformTELEMACinputs_2D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'MARS') && strcmpi(Mode,'3D')
    transformMARSinputs_3D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'MARS') && strcmpi(Mode,'2D')
    transformMARSinputs_2D(conf_name,confOGCM_name)

elseif strcmpi(OGCMmodel_name,'MOHID') && strcmpi(Mode,'3D')
    transformMOHIDinputs_3D(conf_name,confOGCM_name)
    
elseif strcmpi(OGCMmodel_name,'MOHID') && strcmpi(Mode,'2D')
    transformMOHIDinputs_2D(conf_name,confOGCM_name)

end

end