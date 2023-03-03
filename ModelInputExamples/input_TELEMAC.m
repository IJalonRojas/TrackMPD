function conf = input_TELEMAC(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OGCM inputs, folder and files
conf.OGCM.BaseDir = fullfile(pwd,'InputData'); 
conf.OGCM.TELEMACFile = fullfile(conf.OGCM.BaseDir,'hydsed3dNH_ML3_HB5.slf');

% number of points in the new rectangular grid (lon dimension)
conf.OGCM.NumLonGrid = 200; %number of points in the new rectangular grid (lon dimension)

% number of points in the new rectangular grid (lat dimension)
conf.OGCM.NumLatGrid = 200; 

% TELEMAC reference time variable (first time step)
conf.OGCM.t0 = datenum(2018,01,01,00,0,0); 

%%%
