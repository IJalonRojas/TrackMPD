function conf = input_TELEMAC(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OGCM inputs, folder and files
conf.OGCM.BaseDir = fullfile(pwd,'InputData'); 
conf.OGCM.TELEMACFile = fullfile(conf.OGCM.BaseDir,'gir3dNH_sst_frot_3days');

% number of points in the new rectangular grid (lon dimension)
conf.OGCM.NumLonGrid = 500; %number of points in the new rectangular grid (lon dimension)

% number of points in the new rectangular grid (lat dimension)
conf.OGCM.NumLatGrid = 500; 

% TELEMAC reference time variable (first time step)
conf.OGCM.t0 = datenum(2021,06,01,00,0,0); 

%%%
