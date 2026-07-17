function transformInputs(ParameterFile,ModelParameterFile)

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

conf=read_input_conf(ParameterFile);

OGCMmodel_name=conf.OGCM.Model_name;
Mode=conf.Traj.Mode;

if strcmpi(OGCMmodel_name,'POM') && strcmpi(Mode,'3D')
    transformPOMinputs_3D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'POM') && strcmpi(Mode,'2D')
    transformPOMinputs_2D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'FVCOM') && strcmpi(Mode,'3D')
    transformFVCOMinputs_3D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'FVCOM') && strcmpi(Mode,'2D')
    transformFVCOMinputs_2D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'TELEMAC') && strcmpi(Mode,'3D')
% to manage the two options :
    if isfield(conf.OGCM, 'nb_grids') && conf.OGCM.nb_grids == 2
        disp('Starting TELEMAC 3D transformation (Nested 2-Grid Mode)...');
        transformTELEMACinputs_3D_2Grids(conf,ModelParameterFile);
    else
        disp('Starting TELEMAC 3D transformation (Standard 1-Grid Mode)...');
        transformTELEMACinputs_3D_1Grid(conf,ModelParameterFile);
    end
    
elseif strcmpi(OGCMmodel_name,'TELEMAC') && strcmpi(Mode,'2D')
    transformTELEMACinputs_2D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'MARS') && strcmpi(Mode,'3D')
    transformMARSinputs_3D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'MARS') && strcmpi(Mode,'2D')
    transformMARSinputs_2D(conf,ModelParameterFile)

elseif strcmpi(OGCMmodel_name,'MARC')
    transformMARCinputs_2D(conf,ModelParameterFile)

elseif strcmpi(OGCMmodel_name,'MOHID') && strcmpi(Mode,'3D')
    transformMOHIDinputs_3D(conf,ModelParameterFile)
    
elseif strcmpi(OGCMmodel_name,'MOHID') && strcmpi(Mode,'2D')
    transformMOHIDinputs_2D(conf,ModelParameterFile)

end

end
