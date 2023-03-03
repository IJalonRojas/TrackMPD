function conf = input_POM(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%% Particular inputs for POM (SARCCM version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OGCM inputs, folder and files
conf.OGCM.BaseDir = fullfile(pwd,'InputData'); 
conf.OGCM.POM_Prefix = 'u_h';
conf.OGCM.POM_Suffix = '.nc';
    
% POM reference time variable, date of the first time step: datenum(y,m,d,H,M,S)
conf.OGCM.t0 = datenum(1998,06,26,00,0,0); 
   
% Depth for land
conf.OGCM.Hmin = 1; 
       

%%%
