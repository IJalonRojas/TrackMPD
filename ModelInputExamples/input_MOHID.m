function conf = input_MOHID(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Particular inputs for MOHID

% OGCM inputs, folder and files
conf.OGCM.BaseDir = fullfile(pwd,'InputData'); 
conf.OGCM.MOHID_Prefix = 'ext_MOHID_Hydrodynamic_Arousa_';
conf.OGCM.MOHID_Suffix = '.nc';

% Third dimension configuration (only for 3D version)
conf.OGCM.SigmaFile = fullfile(conf.OGCM.BaseDir,'sigma_layers.dat'); 

% Beginning of MOHID time
conf.OGCM.t0 = datenum(2004,01,01,00,00,00); 

% Number of time steps per file    
conf.OGCM.step = 24;  
%%%
    
end